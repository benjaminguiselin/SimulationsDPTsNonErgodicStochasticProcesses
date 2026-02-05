This directory contains everything to create Fig 20.

###################################### 
Simulation codes:
./harmonic_bias_multiple_walkers_parallel.f90
./analysis_harmonic_bias.f90

The code ./harmonic_bias_multiple_walkers_parallel.f90 implements the umbrella sampling procedure discussed in the paper and measures the histogram of the net displacement of all walkers in the presence a harmonic bias (see the documentation of the subroutine INPUT within the code for usage).
The code ./analysis_harmonic_bias.f90 computes the large deviation rate function for the net displacement of all walkers from the different histogram files using the multi-histogram reweighting method (see the documentation of the subroutine INPUT within the code for usage).

###################################### 
Underlying data:
./Data/Batch1/global_params.txt
./Data/Batch1/list_of_simulations.txt
./Data/Batch1/histo*.txt
./Data/Batch1/ldf_3walker.txt
../Data/Batch2/global_params.txt
./Data/Batch2/list_of_simulations.txt
./Data/Batch2/histo*.txt
./Data/Batch2/ldf_3walker.txt
./Data/Batch3/global_params.txt
./Data/Batch3/list_of_simulations.txt
./Data/Batch3/histo*.txt
./Data/Batch3/ldf_3walker.txt

The data are sorted in three folders corresponding to the results for three sets of umbrella sampling potentiels.

The files ./Data/Batch*/global_params.txt show the file of parameters for all umbrella simulations.
The files ./Data/Batch*/list_of_simulations.txt list all umbrella simulations with the following format:
Nber of simulations
Centre of umbrella well         Name of the file where the histogram of end position values are stored

The files ./Data/Batch*/histo*.txt represent the histograms for each umbrella simulation: the 1st row has the centres of histogram bins, the 2nd row has the number of counts.
The files ./Data/Batch*/df_*walker.txt represent the computed large deviation rate function from umbrella simulations.

###################################### 
Figure generation:
./plot_ldf_multiple.gp
./stitch_ldf.py

The script ./stitch_ldf.py collapses the three ldfs computed from the three different batches of simulations by shifting them vertically to generate the file ./Data/ldf_3walker_full.txt.
The gnuplot script ./plot_ldf_multiple.gp generates Fig. 20 from the file ./Data/ldf_3walker_full.txt.