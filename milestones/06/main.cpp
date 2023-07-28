//
// Created by ianni on 6/7/23.
//
#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "lj_neighbors_list.h"
#include "verlet.h"

int main(int argc, char *argv[]) {

    const double epsilon = 1;
    const double sigma = 1;
    const double mass = 1;

    const double neighbour_cutoff = 5;
    const double lj_cutoff = 4.5;
    const double lj_cutoff_engergy = 4 * epsilon * (std::pow(sigma / lj_cutoff, 12) - pow(sigma / lj_cutoff, 6));

    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms = Atoms(names, positions, velocities);
    NeighborList neighbor_list(neighbour_cutoff);
    neighbor_list.update(atoms);


    const double timestep = 0.0003 * std::sqrt(mass * sigma * sigma / epsilon);
    const double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);
    const int num_it = total_time / timestep;
    double e_pot = 0;
    e_pot = lj_neighbors_list(atoms, neighbor_list, epsilon, sigma, lj_cutoff, lj_cutoff_engergy);
    std::cout << "Simulating for at total time of: "<< total_time << "With at timestep of: " << timestep << "\n";
    std::cout << "Resulting in a total of "<< num_it << " itterations.\n";
    std::ofstream traj("milestones/06/traj.xyz");

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
        e_pot = lj_neighbors_list(atoms, neighbor_list, epsilon, sigma, lj_cutoff, lj_cutoff_engergy);
        verlet_step2(atoms, timestep);

        if(i % 1000 == 0){
            std::cout << "Itteration: " << i << "\n";
            write_xyz(traj, atoms);
        }
    }
    traj.close();
}