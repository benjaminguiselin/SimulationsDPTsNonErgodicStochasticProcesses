MODULE kernels

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT64_T, C_DOUBLE

  IMPLICIT NONE
  
  !! PARAMETERS FOR THE PRNG (XOSHIRO256++)
  INTEGER, PARAMETER        :: bits = 64
  INTEGER, PARAMETER        :: state_size = 4
  REAL(C_DOUBLE), PARAMETER :: int64_to_double = 1_C_DOUBLE / REAL((2_C_INT64_T) ** 53, KIND = C_DOUBLE) 

CONTAINS

  !! Initialize the state of the PRNG.
  SUBROUTINE init_state(VectorState)
    INTEGER(C_INT64_T)              :: shift = INT(Z'9E3779B97F4A7C15', KIND = C_INT64_T)
    INTEGER(C_INT64_T)              :: mult1 = INT(Z'BF58476D1CE4E5B9', KIND = C_INT64_T)
    INTEGER(C_INT64_T)              :: mult2 = INT(Z'94D049BB133111EB', KIND = C_INT64_T)
    INTEGER(C_INT64_T)              :: tmp
    INTEGER                         :: i, clock, pid
    INTEGER(C_INT64_T), ALLOCATABLE :: VectorState(:)
    IF (ALLOCATED(VectorState)) DEALLOCATE(VectorState)
    ALLOCATE(VectorState(state_size))
    CALL SYSTEM_CLOCK(COUNT = clock)
    pid = GETPID()
    tmp = ISHFT(INT(clock, KIND = C_INT64_T), 32) + INT(pid, KIND = C_INT64_T)
    DO i = 1, state_size
       tmp = tmp + shift
       tmp = IEOR(tmp, ISHFT(tmp, -30)) * mult1
       tmp = IEOR(tmp, ISHFT(tmp, -27)) * mult2
       tmp = IEOR(tmp, ISHFT(tmp, -31))
       VectorState(i) = tmp
    END DO
    CALL scramble_seed(VectorState)
  END SUBROUTINE init_state
  
  !! Implement a left rotation by k bits for the PRNG.
  FUNCTION rotl(word, k) RESULT(rotation_left)
    INTEGER(C_INT64_T), INTENT(IN) :: word
    INTEGER, INTENT(IN)             :: k
    INTEGER(C_INT64_T)             :: rotation_left
    rotation_left = IOR(ISHFT(word, k), ISHFT(word, k - bits))
  END FUNCTION rotl
  
  !! Update the state of the PRNG.
  SUBROUTINE update_state(VectorState)
    INTEGER(C_INT64_T), INTENT(INOUT) :: VectorState(:)
    INTEGER(C_INT64_T)                :: t
    t = ISHFT(VectorState(2), 17)
    VectorState(3) = IEOR(VectorState(3), VectorState(1))
    VectorState(4) = IEOR(VectorState(4), VectorState(2))
    VectorState(2) = IEOR(VectorState(2), VectorState(3))
    VectorState(1) = IEOR(VectorState(1), VectorState(4))
    VectorState(3) = IEOR(VectorState(3), t)
    VectorState(4) = rotl(VectorState(4), 45)
  END SUBROUTINE update_state
    
  !! Return a random 64-bit integer.
  FUNCTION random_int64(VectorState) RESULT(rand_int)
    INTEGER(C_INT64_T), INTENT(IN) :: VectorState(:)
    INTEGER(C_INT64_T)             :: rand_int
    rand_int = rotl(VectorState(2) * 5, 7) * 9
  END FUNCTION random_int64
  
  !! Return a real number distributed uniformly between 0 (inclusive) and 1 (exclusive).
  FUNCTION random(VectorState) RESULT(rand_nber)
    INTEGER(C_INT64_T), INTENT(INOUT) :: VectorState(:)
    INTEGER(C_INT64_T)                :: rand_int
    REAL(C_DOUBLE)                    :: rand_nber
    rand_int = random_int64(VectorState)
    CALL update_state(VectorState)
    rand_nber = REAL(ISHFT(rand_int, -11), KIND = C_DOUBLE) * int64_to_double
  END FUNCTION random

  !! Scrambe the state of the PRNG in order to avoid poor seeds.
  SUBROUTINE scramble_seed(VectorState)
    INTEGER(C_INT64_T), INTENT(INOUT) :: VectorState(:)
    INTEGER(C_INT64_T)                :: xorkeys(4)
    INTEGER                           :: i
    xorkeys(1) = INT(Z'BD0C5B6E50C2DF49', KIND = C_INT64_T)
    xorkeys(2) = INT(Z'D46061CD46E1DF38', KIND = C_INT64_T)
    xorkeys(3) = INT(Z'BB4F4D4ED6103544', KIND = C_INT64_T)
    xorkeys(4) = INT(Z'114A583D0756AD39', KIND = C_INT64_T)
    DO i = 1, state_size
       VectorState(i) = IEOR(VectorState(i), xorkeys(i))
    END DO
  END SUBROUTINE scramble_seed

  !! Compute the consumed resources for the N walkers up to time T.
  !! eps is the set of all random numbers needed to do the simulation.
  !! For k = 1, N, eps((k-1)*nsteps +1: k*nsteps) follow a uniform distribution between 0 and 1.
  FUNCTION compute_consumed_resources(eps, nwalker, nsteps, dt, T) RESULT(Qcons)
    INTEGER, INTENT(IN)        :: nwalker, nsteps
    REAL(C_DOUBLE), INTENT(IN) :: eps(:), dt, T
    REAL(C_DOUBLE)             :: Qcons, death_rate
    INTEGER                    :: nalive, i, k, index_before
    LOGICAL                    :: alive(nwalker)
    alive = .TRUE.
    index_before = 0
    nalive = nwalker
    death_rate = nalive * dt
    Qcons = 0.0_C_DOUBLE
    DO i = 1, nsteps
       IF (nalive < 1) THEN
          EXIT
       END IF    
       DO k = 1, nwalker
          IF ((alive(k)) .AND. (eps(i + nsteps * ( k - 1 )) < death_rate)) THEN
             alive(k) = .FALSE.
             Qcons = Qcons + nalive * ( i - index_before )
             nalive = nalive - 1
             death_rate = nalive * dt
             index_before = i
          END IF
       END DO
    END DO
    Qcons = Qcons + nalive * ( nsteps - index_before )
    Qcons = Qcons * dt / T
  END FUNCTION compute_consumed_resources
    
  !! Create the entire set eps of random numbers.
  !! For k = 1, N, eps((k-1)*nsteps +1: k*nsteps) follow a uniform distribution between 0 and 1.
  SUBROUTINE create_random_numbers(VectorState, nwalker, nsteps, eps)
    INTEGER(C_INT64_T), INTENT(INOUT) :: VectorState(:)
    INTEGER, INTENT(IN)               :: nwalker, nsteps
    REAL(C_DOUBLE), ALLOCATABLE       :: eps(:)
    INTEGER                           :: i
    IF (.NOT.ALLOCATED(eps)) ALLOCATE(eps(nsteps * nwalker))
    DO i = 1, SIZE(eps)
       eps(i) = random(VectorState)
    END DO
  END SUBROUTINE create_random_numbers
  
  !! Update the array eps of random numbers.
  !! For n_update instants, the random numbers of all particles are updated.
  SUBROUTINE update_random_numbers(VectorState, nwalker, nsteps, n_update, eps)
    INTEGER(C_INT64_T), INTENT(INOUT) :: VectorState(:)
    INTEGER, INTENT(IN)               :: n_update, nsteps, nwalker
    REAL(C_DOUBLE), INTENT(INOUT)     :: eps(:)
    INTEGER                           :: chosen(nsteps), counter, trial_ind, k
    counter = 0
    chosen = 0
    DO WHILE (counter < n_update)
       trial_ind = 1 + INT(nsteps * random(VectorState))
       IF (chosen(trial_ind) == 0) THEN
          chosen(trial_ind) = 1
          counter = counter + 1
          DO k = 1, nwalker
             eps(trial_ind + ( k - 1 ) * nsteps) = random(VectorState)
          END DO
       END IF
    END DO
  END SUBROUTINE update_random_numbers

  !! Read simulation parameters from input files provided on the command line.
  SUBROUTINE input(input_cnf, file_output, file_histo, file_final_state, dt, T, nwalker, &
                   nsteps, ntrials, nsamples, n_update, kappa, x0, xmin, xmax, dx, in_SS)
    CHARACTER(256), INTENT(OUT)     :: input_cnf, file_output, file_final_state, file_histo
    INTEGER, INTENT(OUT)            :: nwalker, nsteps, nsamples, n_update
    INTEGER(C_INT64_T), INTENT(OUT) :: ntrials
    LOGICAL, INTENT(OUT)            :: in_SS
    REAL(C_DOUBLE), INTENT(OUT)     :: dt, T
    REAL(C_DOUBLE), INTENT(OUT)     :: kappa, x0
    REAL(C_DOUBLE), INTENT(OUT)     :: xmin, xmax, dx
    INTEGER                         :: UNIT = 99, in_SS_int
    CHARACTER(256)                  :: global_param_file, local_param_file, param, value_param
    IF (COMMAND_ARGUMENT_COUNT() /= 2) THEN
       WRITE(*, "(a)") "Usage: ./a.out <global_parameter_file> <local_parameter_file>"
       WRITE(*, "(a)") "where <global_parameter_file> is a file with the simulation parameters"
       WRITE(*, "(a)") "for all simulations, and <local_parameter_file> a file with the"
       WRITE(*, "(a)") "simulation parameters particular for this simulation."
       STOP
    END IF
    CALL GET_COMMAND_ARGUMENT(1, VALUE = global_param_file)
    CALL GET_COMMAND_ARGUMENT(2, VALUE = local_param_file)
    input_cnf = ''
    xmin = 1.
    xmax = - 1.
    dx = -1.
    in_SS_int = -1
    nwalker = 0
    OPEN(UNIT = UNIT, FILE = global_param_file, STATUS = "OLD")
    DO
       READ(UNIT, *, END = 10) param, value_param
       IF (TRIM(param) == "dt") READ(value_param(1:LEN_TRIM(value_param)),*) dt !! time step for the resolution of Newton's equations
       IF (TRIM(param) == "nsteps") READ(value_param(1:LEN_TRIM(value_param)),*) nsteps !! total number of steps for one trajectory
       IF (TRIM(param) == "ntrials") READ(value_param(1:LEN_TRIM(value_param)),*) ntrials !! total number of simulated trajectories
       IF (TRIM(param) == "n_update") READ(value_param(1:LEN_TRIM(value_param)),*) n_update !! number of random numbers to change between two states of the MC
       IF (TRIM(param) == "nsamples") READ(value_param(1:LEN_TRIM(value_param)),*) nsamples !! number of samples to store if not in the steady state
       IF (TRIM(param) == "xmin") READ(value_param(1:LEN_TRIM(value_param)),*) xmin !! minimal bin value for the histogram
       IF (TRIM(param) == "xmax") READ(value_param(1:LEN_TRIM(value_param)),*) xmax !! maximal bin value for the histogram
       IF (TRIM(param) == "dx") READ(value_param(1:LEN_TRIM(value_param)),*) dx !! bin width for the histogram
       IF (TRIM(param) == "kappa") READ(value_param(1:LEN_TRIM(value_param)),*) kappa !! spring constant for the harmonic potential
       IF (TRIM(param) == "nwalker") READ(value_param(1:LEN_TRIM(value_param)),*) nwalker !! number of walkers at the beginning of the simulation
    END DO
10  CONTINUE
    CLOSE(UNIT = UNIT)
    OPEN(UNIT = UNIT, FILE = local_param_file, STATUS = "OLD")
    DO
       READ(UNIT, *, END = 17) param, value_param
       IF (TRIM(param) == "ntrials") READ(value_param(1:LEN_TRIM(value_param)),*) ntrials !! total number of simulated trajectories
       IF (TRIM(param) == "n_update") READ(value_param(1:LEN_TRIM(value_param)),*) n_update !! number of random numbers to change between two states of the MC
       IF (TRIM(param) == "file_output") file_output = value_param !! name of file where to store the samples
       IF (TRIM(param) == "file_histo") file_histo = value_param !! name of file where to store the histogram
       IF (TRIM(param) == "file_input") input_cnf = value_param !! name of file where to find the initial state (if provided)
       IF (TRIM(param) == "file_final_state") file_final_state = value_param !! name of file where to store the final state
       IF (TRIM(param) == "x0") READ(value_param(1:LEN_TRIM(value_param)),*) x0 !! center of the harmonic bias
       IF (TRIM(param) == "in_SS") READ(value_param(1:LEN_TRIM(value_param)),*) in_SS_int !! 1 if in_SS (and compute histograms), 0 otherwise    
    END DO
17  CONTINUE
    CLOSE(UNIT = UNIT)
    T = dt * nsteps
    IF (in_SS_int .EQ. 1) THEN
       in_SS = .TRUE.
    ELSE IF (in_SS_int .EQ. 0) THEN
       in_SS = .FALSE.
    ELSE
       WRITE(*, "(a)") "in_SS should be 0 or 1."
       STOP
    END IF
    IF ( (in_SS) .AND. ( (dx < 0) .OR. (xmin > xmax)) ) THEN
       WRITE(*, "(a)") "Range and bin width must be provided to compute the histogram."
       STOP
    END IF
    IF (nwalker < 1) THEN
       WRITE(*, "(a)") "The number of walkers must be strictly positive."
       STOP
    END IF
  END SUBROUTINE input
  
  !! Read the list of all random numbers from a file.
  SUBROUTINE read_init_state(UNIT, nwalker, nsteps, eps)
    INTEGER, INTENT(IN)         :: UNIT, nwalker, nsteps
    REAL(C_DOUBLE), ALLOCATABLE :: eps(:)
    IF (.NOT.ALLOCATED(eps)) ALLOCATE(eps(nwalker * nsteps))
    READ(UNIT, *) eps(:)
  END SUBROUTINE read_init_state
  
  !! Write the list of all random numbers in a file.
  SUBROUTINE write_state(UNIT, eps)
    INTEGER, INTENT(IN)        :: UNIT
    REAL(C_DOUBLE), INTENT(IN) :: eps(:)
    WRITE(UNIT, *) eps(:)
  END SUBROUTINE write_state
  
  !! Report the time series of samples.
  SUBROUTINE report_samples(UNIT, list_x)
    INTEGER, INTENT(IN)        :: UNIT
    REAL(C_DOUBLE), INTENT(IN) :: list_x(:)
    INTEGER                    :: i
    DO i = 1, SIZE(list_x)
       WRITE(UNIT, "(ES20.8E3)") list_x(i)
    END DO
  END SUBROUTINE report_samples
  
  !! Add the value x to the histogram.
  SUBROUTINE add_at_histogram(x, bin_edges, histo)
    INTEGER(C_INT64_T), INTENT(INOUT) :: histo(:)
    REAL(C_DOUBLE), INTENT(IN)        :: x, bin_edges(:)
    INTEGER                           :: i
    i = MINLOC(bin_edges, MASK = bin_edges .GE. x, DIM = 1)
    i = i - 1
    IF (i .GE. 1) THEN
        histo(i) = histo(i) + 1
    END IF
  END SUBROUTINE add_at_histogram
  
  !! Report the histogram of values.
  SUBROUTINE report_histo(UNIT, bin_edges, histo)
    INTEGER, INTENT(IN)            :: UNIT
    INTEGER(C_INT64_T), INTENT(IN) :: histo(:)
    REAL(C_DOUBLE), INTENT(IN)     :: bin_edges(:)
    INTEGER                        :: i
    DO i = 1, SIZE(histo)
       WRITE(UNIT, "(F8.4, I20)") 0.5 * ( bin_edges(i) + bin_edges(i + 1) ), histo(i)
    END DO
  END SUBROUTINE report_histo

END MODULE kernels

PROGRAM main
  
  USE kernels
  
  IMPLICIT NONE
  
  REAL(C_DOUBLE), ALLOCATABLE     :: eps(:), eps_trial(:)
  REAL(C_DOUBLE), ALLOCATABLE     :: list_x(:), bin_edges(:)
  INTEGER(C_INT64_T), ALLOCATABLE :: VectorState(:)
  INTEGER(C_INT64_T), ALLOCATABLE :: histo(:)
  REAL(C_DOUBLE)                  :: dt, T
  REAL(C_DOUBLE)                  :: kappa, x0
  REAL(C_DOUBLE)                  :: xmin, xmax, dx
  REAL(C_DOUBLE)                  :: trial_x, bias, x
  REAL(C_DOUBLE)                  :: start_time, save_time
  INTEGER(C_INT64_T)              :: ntrials, i, freq_save
  INTEGER                         :: nwalker, nsteps, n_update
  INTEGER                         :: nsamples, nbin, j
  INTEGER                         :: idx
  INTEGER                         :: unit_init = 101, unit_output = 102
  INTEGER                         :: unit_final_state = 103, unit_histo = 104
  CHARACTER(256)                  :: file_output, file_final_state, file_init, file_histo
  LOGICAL                         :: acc_move
  LOGICAL                         :: in_SS
  
  !! Get global parameters.
  CALL input(file_init, file_output, file_histo, file_final_state, dt, T, nwalker, nsteps, &
             ntrials, nsamples, n_update, kappa, x0, xmin, xmax, dx, in_SS)  
  freq_save = MAX(ntrials / INT(nsamples, KIND = C_INT64_T), 1_C_INT64_T)
               
  !! Initialize the random numbers.
  CALL init_state(VectorState)
  
  !! Create the first array of random numbers or read it from the initial state.
  IF (LEN_TRIM(file_init) > 0) THEN
     OPEN(UNIT = unit_init, FILE = file_init, STATUS = "UNKNOWN")
     CALL read_init_state(unit_init, nwalker, nsteps, eps)
     CLOSE(UNIT = unit_init)
  ELSE
     CALL create_random_numbers(VectorState, nwalker, nsteps, eps)
  END IF
  
  !! Create the array of samples.
  ALLOCATE(list_x(nsamples))
  idx = 0
  
  !! Create the histogram and the bin edges if in the steady state.
  IF (in_SS) THEN
     nbin = INT(( xmax - xmin ) / dx) + 1
     ALLOCATE(bin_edges(nbin + 1), histo(nbin))
     histo = 0
     DO j = 1, nbin + 1
        bin_edges(j) = xmin + ( j - 1 ) * dx
     END DO
  END IF
  
  !! Compute the different trajectories.
  x = compute_consumed_resources(eps, nwalker, nsteps, dt, T)
  ALLOCATE(eps_trial(SIZE(eps)))
  CALL CPU_TIME(start_time)
  DO i = 1_C_INT64_T, ntrials
     eps_trial = eps
     CALL update_random_numbers(VectorState, nwalker, nsteps, n_update, eps_trial)
     trial_x = compute_consumed_resources(eps_trial, nwalker, nsteps, dt, T)
     bias = ( trial_x - x0 ) ** 2 - ( x - x0 ) ** 2
     IF (bias < 0.) THEN
        acc_move = .TRUE.
     ELSE IF (random(VectorState) < EXP(- kappa * bias)) THEN
        acc_move = .TRUE.
     ELSE
        acc_move = .FALSE.
     END IF
     IF (acc_move) THEN
        x = trial_x
        eps = eps_trial
     END IF
     IF (in_SS) THEN
        CALL add_at_histogram(x, bin_edges, histo)
     END IF
     IF (MODULO(i, freq_save) < 1) THEN
        idx = idx + 1
        list_x(idx) = x
     END IF
     CALL CPU_TIME(save_time)
     IF (save_time - start_time > 28800) THEN
        OPEN(UNIT = unit_final_state, FILE = file_final_state, STATUS = "UNKNOWN")
        CALL write_state(unit_final_state, eps)
        CLOSE(UNIT = unit_final_state)
        OPEN(UNIT = unit_output, FILE = file_output, STATUS = "UNKNOWN")
        CALL report_samples(unit_output, list_x(1:idx))
        CLOSE(UNIT = unit_output)
        IF (in_SS) THEN
           OPEN(UNIT = unit_histo, FILE = file_histo)
           CALL report_histo(unit_histo, bin_edges, histo)
           CLOSE(UNIT = unit_histo)
        END IF
        start_time = save_time
     END IF
  END DO
  
  !! Store the final list of random numbers.
  OPEN(UNIT = unit_final_state, FILE = file_final_state, STATUS = "UNKNOWN")
  CALL write_state(unit_final_state, eps)
  CLOSE(UNIT = unit_final_state)
  
  !! Store the time series of samples.
  OPEN(UNIT = unit_output, FILE = file_output, STATUS = "UNKNOWN")
  CALL report_samples(unit_output, list_x(1:idx))
  CLOSE(UNIT = unit_output)
  
  !! Report the histogram.
  IF (in_SS) THEN
     OPEN(UNIT = unit_histo, FILE = file_histo)
     CALL report_histo(unit_histo, bin_edges, histo)
     CLOSE(UNIT = unit_histo)
  END IF
  
END PROGRAM main