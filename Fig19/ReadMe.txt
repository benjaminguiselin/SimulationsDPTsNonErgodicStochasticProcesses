This directory contains everything to create Fig 19.

###################################### 
Simulation codes:
./harmonic_bias_parallel.f90
./analysis_harmonic_bias.f90

The code ./harmonic_bias_parallel.f90 implements the umbrella sampling procedure discussed in the paper and measures the histogram of the net resources consumed by the living walkers in the presence a harmonic bias (see the documentation of the subroutine INPUT within the code for usage).
The code ./analysis_harmonic_bias.f90 computes the large deviation rate function for the net resources consumed by the living walkers from the different histogram files using the multi-histogram reweighting method (see the documentation of the subroutine INPUT within the code for usage).

###################################### 
Underlying data:
./Data/One_walker/global_params_1walker.txt
./Data/One_walker/list_of_simulations_1walker.txt
./Data/One_walker/histo*.txt
./Data/One_walker/ldf_1walker.txt
./Data/Two_walkers/global_params_2walker.txt
./Data/Two_walkers/list_of_simulations_2walker.txt
./Data/Two_walkers/histo*.txt
./Data/Two_walkers/ldf_2walker.txt
./Data/Three_walkers/global_params_3walker.txt
./Data/Three_walkers/list_of_simulations_3walker.txt
./Data/Three_walkers/histo*.txt
./Data/Three_walkers/ldf_3walker.txt

The data are sorted in three folders corresponding to the results for one walker (./Data/One_walker/), two walkers (./Data/Two_walkers/) and three walkers (./Data/Three_walkers/).

The files ./Data/*_walker*/global_params_*walker.txt show the file of parameters for all umbrella simulations.
The files ./Data/*_walker*/list_of_simulations_*walker.txt list all umbrella simulations with the following format:
Nber of simulations
Centre of umbrella well         Name of the file where the histogram of end position values are stored

The files ./Data/*_walker*/histo*.txt represent the histograms for each umbrella simulation: the 1st row has the centres of histogram bins, the 2nd row has the number of counts.
The files ./Data/*_walker*/ldf_*walker.txt represent the computed large deviation rate function from umbrella simulations.

###################################### 
Figure generation:
./plot_ldf_epidemiology.gp

The gnuplot script ./plot_ldf_epidemiology.gp generates Fig. 19 from the files ./Data/One_walker/ldf_1walker.txt, ./Data/Two_walkers/ldf_2walker.txt and ./Data/Three_walkers/ldf_3walker.txt.
