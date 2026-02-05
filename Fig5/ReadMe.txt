This directory contains everything to create Fig 5.

###################################### 
Simulation codes:
./harmonic_bias_residence_time_parallel.f90
./analysis_harmonic_bias.f90

The code ./harmonic_bias_residence_time_parallel.f90 implements the umbrella sampling procedure discussed in the paper and measures the histogram of the residence time of the walker in the presence a harmonic bias (see the documentation of the subroutine INPUT within the code for usage).
The code ./analysis_harmonic_bias.f90 computes the large deviation rate function for the residence time of the walker from the different histogram files using the multi-histogram reweighting method (see the documentation of the subroutine INPUT within the code for usage).

###################################### 
Underlying data:
./Data/global_params.txt
./Data/list_of_simulations.txt
./Data/histo*.txt
./Data/ldf_residence_time.txt

The file ./Data/global_params.txt shows the file of parameters for all umbrella simulations.
The file ./Data/list_of_simulations.txt lists all umbrella simulations with the following format:
Nber of simulations
Centre of umbrella well         Name of the file where the histogram of end position values are stored

The files ./Data/histo*.txt represent the histograms for each umbrella simulation: the 1st row has the centres of histogram bins, the 2nd row has the number of counts.
The file ./Data/ldf_residence_time.txt represents the computed large deviation rate function from umbrella simulations.

###################################### 
Figure generation:
./plot_ldf_residence_time.gp
./compute_ldf_residence_time.py

The script ./compute_ldf_residence_time.py computes the exact analytical solution for the large deviation rate function (./Analytics/ldf_residence_time_theory.txt).
The gnuplot script ./plot_ldf_residence_time.gp generates Fig. 5 from the files ./Data/ldf_residence_time.txt and ./Analytics/ldf_residence_time_theory.txt.