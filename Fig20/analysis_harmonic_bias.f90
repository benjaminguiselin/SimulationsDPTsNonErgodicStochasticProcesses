MODULE kernels_analysis

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT64_T, C_DOUBLE, C_LONG_DOUBLE
  
  IMPLICIT NONE
 
  REAL(C_LONG_DOUBLE), PARAMETER :: delta = 1E-7 !! Criterion to stop the iteration for computing the partition functions

CONTAINS

  !! Read simulation parameters from input files provided on the command line.
  SUBROUTINE input(file_histo, file_ldf, nsimu, kappa, x0, xmin, xmax, dx, dx_new, T)
    CHARACTER(256), ALLOCATABLE :: file_histo(:)
    CHARACTER(256), INTENT(OUT) :: file_ldf
    INTEGER, INTENT(OUT)        :: nsimu
    REAL(C_DOUBLE), ALLOCATABLE :: x0(:)
    REAL(C_DOUBLE), INTENT(OUT) :: kappa, T
    REAL(C_DOUBLE), INTENT(OUT) :: xmin, xmax, dx, dx_new
    INTEGER                     :: UNIT = 99, i, nsteps
    REAL(C_DOUBLE)              :: dt
    CHARACTER(256)              :: global_param_file, list_of_simulations, param, value_param
    IF (COMMAND_ARGUMENT_COUNT() /= 2) THEN
       WRITE(*, "(a)") "Usage: ./analysis.out <global_parameter_file> <list_of_simulations>"
       WRITE(*, "(a)") "where <global_parameter_file> is a file with the simulation parameters"
       WRITE(*, "(a)") "for all simulations, and <list_of_simulations> a file with the"
       WRITE(*, "(a)") "list of all umbrella simulations."
       STOP
    END IF
    CALL GET_COMMAND_ARGUMENT(1, VALUE = global_param_file)
    CALL GET_COMMAND_ARGUMENT(2, VALUE = list_of_simulations)
    OPEN(UNIT = UNIT, FILE = global_param_file, STATUS = "OLD")
    dx_new = 0.0_C_DOUBLE
    DO
       READ(UNIT, *, END = 10) param, value_param
       IF (TRIM(param) == "xmin") READ(value_param(1:LEN_TRIM(value_param)),*) xmin !! minimal bin value for the histogram
       IF (TRIM(param) == "xmax") READ(value_param(1:LEN_TRIM(value_param)),*) xmax !! maximal bin value for the histogram
       IF (TRIM(param) == "dx") READ(value_param(1:LEN_TRIM(value_param)),*) dx !! bin width for the histogram
       IF (TRIM(param) == "dx_new") READ(value_param(1:LEN_TRIM(value_param)),*) dx_new !! new bin width for the histogram (if coarsening is necessary)
       IF (TRIM(param) == "kappa") READ(value_param(1:LEN_TRIM(value_param)),*) kappa !! spring constant for the harmonic potential
       IF (TRIM(param) == "file_ldf") file_ldf = value_param !! name of the file where to store the numerical ldf
       IF (TRIM(param) == "dt") READ(value_param(1:LEN_TRIM(value_param)),*) dt !! time step for the resolution of Newton's equations
       IF (TRIM(param) == "nsteps") READ(value_param(1:LEN_TRIM(value_param)),*) nsteps !! total number of steps for one trajectory
    END DO
10  CONTINUE
    CLOSE(UNIT = UNIT)
    T = dt * nsteps
    dx_new = INT(MAX(dx, dx_new) / dx) * dx
    OPEN(UNIT = UNIT, FILE = list_of_simulations, STATUS = "OLD")
    READ(UNIT, *) nsimu
    ALLOCATE(file_histo(nsimu)) !! List of the names of all histogram files
    ALLOCATE(x0(nsimu)) !! List of all centres of the umbrella potentials
    DO i = 1, nsimu
       READ(UNIT, *) x0(i), file_histo(i)  
    END DO
    CLOSE(UNIT = UNIT)
  END SUBROUTINE input
  
  !! Report the numerical ldf.
  SUBROUTINE report_ldf(UNIT, T, bin_centers, pdf)
    INTEGER                         :: UNIT
    REAL(C_DOUBLE), INTENT(IN)      :: bin_centers(:), T
    REAL(C_LONG_DOUBLE), INTENT(IN) :: pdf(:)
    REAL(C_LONG_DOUBLE)             :: max_val
    INTEGER                         :: j
    max_val = LOG(MAXVAL(pdf))
    DO j = 1, SIZE(pdf)
       IF (pdf(j) > 0.0_C_LONG_DOUBLE) THEN
          WRITE(UNIT, "(F8.4, ES20.8E3)") bin_centers(j), ( max_val - LOG(pdf(j)) ) / T
       END IF
    END DO
  END SUBROUTINE report_ldf
  
  !! Read the histogram from file.
  SUBROUTINE read_histo(UNIT, nbin, ncoarsen, histo, nsamples)
    INTEGER, INTENT(IN)             :: UNIT, nbin, ncoarsen
    INTEGER(C_INT64_T), INTENT(OUT) :: histo(:), nsamples
    INTEGER                         :: j, jj
    INTEGER(C_INT64_T)              :: histo_tmp
    REAL(C_DOUBLE)                  :: tmp
    nsamples = 0_C_INT64_T
    DO j = 1, nbin
       histo(j) = 0_C_INT64_T
       DO jj = ( j - 1 ) * ncoarsen + 1, j * ncoarsen
          READ(UNIT, *) tmp, histo_tmp
          histo(j) = histo(j) + histo_tmp
       END DO
       nsamples = nsamples + histo(j)
    END DO
  END SUBROUTINE read_histo
  
END MODULE kernels_analysis

PROGRAM main
  
  USE kernels_analysis
  
  IMPLICIT NONE
  
  REAL(C_DOUBLE), ALLOCATABLE      :: x0(:)
  REAL(C_DOUBLE), ALLOCATABLE      :: bin_edges(:), bin_centers(:)
  REAL(C_LONG_DOUBLE), ALLOCATABLE :: pdf(:), part_func(:), part_func_before(:)
  REAL(C_LONG_DOUBLE), ALLOCATABLE :: log_part_func(:), nsamples_real(:)
  REAL(C_LONG_DOUBLE), ALLOCATABLE :: numerator_WHAM(:), bias(:,:), log_nsamp(:)
  INTEGER(C_INT64_T), ALLOCATABLE  :: histo(:,:), nsamples(:)
  CHARACTER(256), ALLOCATABLE      :: file_histo(:)
  REAL(C_DOUBLE)                   :: kappa, T
  REAL(C_DOUBLE)                   :: xmin, xmax, dx, dx_new
  REAL(C_LONG_DOUBLE)              :: max_arg_exp, log_denominator
  INTEGER                          :: nsimu, nbin
  INTEGER                          :: i, j, imax
  INTEGER                          :: unit_histo = 101, unit_ldf = 102
  CHARACTER(256)                   :: file_ldf
  
  !! Get global parameters.
  CALL input(file_histo, file_ldf, nsimu, kappa, x0, xmin, xmax, dx, dx_new, T)

  !! Create the histogram, the bin edges and the bin centers.
  nbin = INT(( xmax - xmin ) / dx_new) + 1
  IF (xmin + ( INT(( xmax - xmin ) / dx_new) + 1 ) * dx_new > &
  xmin + dx * ( INT(( xmax - xmin ) / dx) + 1 )) THEN
     nbin = nbin - 1
  END IF
  ALLOCATE(bin_edges(nbin + 1), bin_centers(nbin))
  ALLOCATE(histo(nbin, nsimu), nsamples(nsimu))
  DO j = 1, nbin + 1
     bin_edges(j) = xmin + ( j - 1 ) * dx_new
  END DO
  DO j = 1, nbin
     bin_centers(j) = 0.5 * ( bin_edges(j) + bin_edges(j + 1) )
  END DO
  
  !! Read the histograms from the files.
  DO i = 1, nsimu
     OPEN(UNIT = unit_histo, FILE = file_histo(i), STATUS = "OLD")
     CALL read_histo(unit_histo, nbin, INT(dx_new / dx), histo(:,i), nsamples(i))
     CLOSE(UNIT = unit_histo)
  END DO
  
  !! Compute the partition functions.
  ALLOCATE(part_func(nsimu), part_func_before(nsimu))
  part_func_before = 2.0_C_LONG_DOUBLE
  part_func = 1.0_C_LONG_DOUBLE
  ALLOCATE(nsamples_real(nsimu), log_part_func(nsimu), log_nsamp(nsimu))
  nsamples_real = REAL(nsamples, KIND = C_LONG_DOUBLE)
  log_nsamp = LOG(nsamples_real)
  ALLOCATE(numerator_WHAM(nbin), bias(nbin, nsimu))
  DO j = 1, nbin
     numerator_WHAM(j) =  REAL(SUM(histo(j,:)), KIND = C_LONG_DOUBLE)
     DO i = 1, nsimu
        bias(j,i) = REAL(kappa * ( bin_centers(j) - x0(i) ) ** 2, KIND = C_LONG_DOUBLE)
     END DO
  END DO
  DO WHILE (MAXVAL(ABS(part_func_before / part_func - 1.0_C_LONG_DOUBLE)) > delta)
     part_func_before = part_func
     log_part_func = log_nsamp - LOG(part_func)
     part_func = 0.0_C_LONG_DOUBLE
     DO j = 1, nbin
        DO i = 1, nsimu
           imax = MAXLOC(log_part_func - bias(j, :), DIM = 1)
           max_arg_exp = log_part_func(imax) - bias(j,imax)
           log_denominator = max_arg_exp + bias(j,i) + LOG(SUM(EXP(log_part_func - bias(j,:) - max_arg_exp)))
           part_func(i) = part_func(i) + numerator_WHAM(j) * EXP(- log_denominator)
        END DO
     END DO
     !PRINT*, part_func
     part_func = part_func / SQRT(MINVAL(part_func) * MAXVAL(part_func))
     PRINT*, MAXVAL(ABS(part_func_before / part_func - 1.))
  END DO
  
  !! Compute the best estimate of the ldf.
  ALLOCATE(pdf(nbin))
  DO j = 1, nbin
     pdf(j) = numerator_WHAM(j) / SUM(nsamples_real * EXP(- bias(j,:)) / part_func)
  END DO
  OPEN(UNIT = unit_ldf, FILE = file_ldf, STATUS = "UNKNOWN")
  CALL report_ldf(unit_ldf, T, bin_centers, pdf)
  CLOSE(UNIT = unit_ldf)
    
END PROGRAM main
