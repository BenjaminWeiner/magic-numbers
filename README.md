# magic-numbers
Code for 3D lattice simulations in the article, "Rigidity enhances a magic-number effect in polymer phase separation" (https://arxiv.org/abs/1901.09352)

******************************************
Simulation executables:

Both the Windows and Linux executables accept the following arguments:
test name (labels output),
lattice specification file (size, connectivity), 
polymer specification file (length and number of flexible polymers, number of rigid cubes),
temperature schedule for simulated annealing,
Number of Monte Carlo steps (one Monte Carlo step is x moves, where x=# of monomers in simulation),
Strength of nonspecific interactions, in k_BT,

For example, the following commands would run the Windows and Linux versions, respectively:
magicnumbers_3D.exe testName lattice_cubic.txt polySpecs_lP7_c0.09 annealing_schedule_4 40000 -0.1
./magicnumbers_3D_linux testName lattice_cubic.txt polySpecs_lP7_c0.09 annealing_schedule_4 40000 -0.1

The output will be the polyOutput_testName.dat, which has the position of every monomer at the conclusion of each Monte Carlo step. All the flexible polymers are listed before all the rigid blocks.
This output can be analyzed using Python script "magicnumbers_clustersize.py," described below. 
*************************************************************
Analysis script: magicnumbers_clustersize.py

This script extracts the sizes of connected clusters of polymers from the simulation output. It saves the list of cluster sizes (in monomers) at each time step.
The only argument is "testName," which you used to label the output from the executable. Note that this code requires the "params" record file produced by the simulation as well.
