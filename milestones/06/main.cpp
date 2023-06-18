//
// Created by ianni on 6/7/23.
//
#include <iostream>
#include "xyz.h"
int main(int argc, char *argv[]) {

    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms = Atoms(names, positions, velocities);

}