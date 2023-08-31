#! /bin/bash


executable1=cmake-build-release/final/simulation_1/simulation_1_EAMCutoff
executable2=cmake-build-release/final/simulation_1/simulation_1_LatticeConstant
executable3=cmake-build-release/final/simulation_1/simulation_1_Timestep

mpirun -n 8 $executable1
mpirun -n 8 $executable2
mpirun -n 8 $executable3