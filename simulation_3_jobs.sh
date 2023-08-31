#! /bin/bash


executable=cmake-build-release/final/simulation_3/simulation_3

# Different initial temperatures
#                       dt cutoff radius strain-rate% max-strain%  initial-temperature
mpirun -n 8 $executable 10 7      14     0.00001      10           000
mpirun -n 8 $executable 10 7      14     0.00001      10           100
mpirun -n 8 $executable 10 7      14     0.00001      10           200
mpirun -n 8 $executable 10 7      14     0.00001      10           300
mpirun -n 8 $executable 10 7      14     0.00001      10           400
mpirun -n 8 $executable 10 7      14     0.00001      10           500
mpirun -n 8 $executable 10 7      14     0.00001      10           600
mpirun -n 8 $executable 10 7      14     0.00001      10           700
mpirun -n 8 $executable 10 7      14     0.00001      10           800

# Different radii
#                       dt cutoff radius strain-rate% max-strain%  initial-temperature
mpirun -n 8 $executable 10 7      10     0.00001      10           300
mpirun -n 8 $executable 10 7      11     0.00001      10           300
mpirun -n 8 $executable 10 7      12     0.00001      10           300
mpirun -n 8 $executable 10 7      13     0.00001      10           300
mpirun -n 8 $executable 10 7      15     0.00001      10           300
mpirun -n 8 $executable 10 7      17     0.00001      10           300
mpirun -n 8 $executable 10 7      19     0.00001      10           300

# Different strain-rate
#                       dt cutoff radius strain-rate% max-strain%  initial-temperature
mpirun -n 8 $executable 10 7      14     0.000005     10           300
mpirun -n 8 $executable 10 7      14     0.000006     10           300
mpirun -n 8 $executable 10 7      14     0.000007     10           300
mpirun -n 8 $executable 10 7      14     0.000008     10           300
mpirun -n 8 $executable 10 7      14     0.000009     10           300
mpirun -n 8 $executable 10 7      14     0.000010     10           300
mpirun -n 8 $executable 10 7      14     0.000100     10           300