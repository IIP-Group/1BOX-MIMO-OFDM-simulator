# Simulator for 1-Bit OFDM (1BOX) detector for massive MU-MIMO-OFDM
(c) 2020 Seyed Hadi Mirfarshbafan
e-mail: sm2675@cornell.edu 


### Important information 

If you are using the simulator (or parts of it) for a publication, then you *must* cite our paper:

S. H. Mirfarshbafan, M. Shabany, S. A. Nezamalhosseini, and C. Studer, "Algorithm and VLSI Design for 1-bit Data Detection in Massive MIMO-OFDM", IEEE Open J. of Circuits and Sys. (OJCAS), to appear in 2020

and clearly mention this in your paper.  In the following, by *paper*, we mean the above paper.

### How to start a simulation:

To regenerate the plots in figure 3 of the paper, simply assign one of the letters 'a', 'b', 'c' or 'd' to the variable *simulation_setup* in the main_BER_sim.m 
and run the script. It will generate all curves corresponding to the chosen letter from figure 3 and plot them all, except for the fixed-point version of 1BOX detector.

In order to simulate with other parameters, define a new set of parameters in par_config.m and assign a value to the variable simulation_setup that corresponds to the
new set of parameters. Additionally, define what combination of channel estimator and data detectors to be simulated, in the par.chest_list and par.det_list cell arrays
in the par_config.m file. The cell array par.quant_list, specifies for each pair of channel estimator/data detector, whether infinite resolution or 1-bit ADCs
are considered.
 	
We highly recommend you to execute the code step-by-step (using MATLAB's debug mode) in order to get a detailed understanding of the simulator. In order to do so, you should change the parfor in the main_BER_sim.m to a simple for (just replace parfor with for) so that debugging is possible.

### Version history
* Version 1 (Sep. 03 2020) - sm2675@cornell.edu  - uploaded the simulation
