#! /bin/bash


executable=cmake-build-release/final/simulation_2/simulation_2

#                       cluster dt  cutoff create-trajectory-file itterations-pre-frame start-T end-T Delta-T tau_relax
mpirun -n 8 $executable 923     10  7      true                   500                   600     800   2       2000
mpirun -n 8 $executable 1415    10  7      true                   500                   600     850   2       2000
mpirun -n 8 $executable 2057    10  7      true                   500                   600     850   2       2000
mpirun -n 8 $executable 2869    10  7      true                   500                   600     850   2       2000
mpirun -n 8 $executable 3871    10  7      true                   500                   600     850   2       2000
mpirun -n 8 $executable 5083    10  7      true                   500                   600     850   2       2000
mpirun -n 8 $executable 6525    10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 8217    10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 10179   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 12431   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 14993   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 17885   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 21127   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 24739   10  7      true                   500                   600     1100  2       2000
mpirun -n 8 $executable 28741   10  7      true                   500                   600     1100  2       2000