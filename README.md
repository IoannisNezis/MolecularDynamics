# Molecular Dynamics Code for IMTEK-HPC-MD class

## Usage
The final simulation code is in the folder `final`.
There you will find 3 directories, one for each simulation described in the report.
To run the simulations build the targets:
- simulation_1_Timestep
- simulation_1_EAMCutoff
- simulation_1_LatticeConstant
- simulation_2
- simulation_3

Remember to set the Working directory and to enable `-DUSE_MPI=true` in the CMake profile `Release`.

Then you can run the simulations.

The simulations described in the report can be reproduced by running the shell scripts:
- simulation_1_jobs.sh
- simulation_2_jobs.sh
- simulation_3_jobs.sh

Each simulation consumes input files from their `input` directory and writes the data and trajectory files to:
- `output/data`
- `output/trajectory`
