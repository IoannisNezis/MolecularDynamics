//
// Created by ianni on 5/17/23.
//
#include <iostream>
#include "lj_direct_summation.h"
#include "xyz.h"
#include "verlet.h"


int main(int argc, char *argv[]) {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms = Atoms(positions, velocities);
    const double epsilon = 1;
    const double sigma = 1;
    const double mass = 1;
    const double timestep = 0.001 * std::sqrt(mass * sigma * sigma / epsilon);
    const double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);
    const int num_it = total_time / timestep;
    double e_pot;
    e_pot = lj_direct_summation(atoms);
    std::cout << "Simulating for at total time of: "<< total_time << "With at timestep of: " << timestep << "\n";
    std::cout << "Resulting in a total of "<< num_it << " itterations.\n";
    std::ofstream traj("traj.xyz");

    for (int i = 0; i < num_it; ++i) {
//        std::cout << "Epot: " << e_pot << "\n";
//        std::cout << "Ekin: " << atoms.kinetic_engergy() << "\n";
//        std::cout << "E   : " << e_pot + atoms.kinetic_engergy() << "\n";
//
//        std::cout << "---------- First position -----------"<< "\n";
//        std::cout << atoms.positions.col(0) << "\n";
//        std::cout << "---------- First velocity -----------"<< "\n";
//        std::cout << atoms.velocities.col(0) << "\n";
//        std::cout << "---------- First force -----------"<< "\n";
//        std::cout << atoms.forces.col(0) << "\n";

        verlet_step1(atoms, timestep);
        e_pot = lj_direct_summation(atoms);
        verlet_step2(atoms, timestep);

        if(i % 1000 == 0){
            std::cout << "Itteration: " << i << "\n";
            write_xyz(traj, atoms);
        }
    }
    traj.close();
}