//
// Created by ianni on 7/5/23.
//
#include <iostream>
#include <iomanip>
#include "xyz.h"
#include "neighbors.h"
#include "lj_neighbors_list.h"
#include "verlet.h"
#include "thermostat.h"

void create_lattice(int nx,int ny,int nz,double lattice_constant){
    std::ofstream file("milestones/05/lattice.xyz");
    file << nx*ny*nz << std::endl << std::endl;
    float dx,dy,dz;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                dx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                dy = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                dz = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                file << std::setw(2) << "H" << " "
                     << std::setw(10) << x*lattice_constant << " " << y*lattice_constant << " " << z*lattice_constant
                     << std::setw(13) << dx << " " << dy << " " << dz
                     << std::endl;
            }
        }
    }
    file.close();
}

int main(int argc, char *argv[]) {
    const double epsilon = 1;
    const double sigma = 1;
    const double mass = 1;

    const double neighbour_cutoff = 510.5;
    const double lj_cutoff = 510;
    const double lj_cutoff_engergy = 4 * epsilon * (std::pow(sigma / lj_cutoff, 12) - pow(sigma / lj_cutoff, 6));

    create_lattice(5,5,5,2* sigma);

    auto [names, positions, velocities]{read_xyz_with_velocities("milestones/05/lattice.xyz")};
    Atoms atoms = Atoms(names, positions, velocities);
    NeighborList neighbor_list(neighbour_cutoff);
    neighbor_list.update(atoms);
    atoms.forces.setZero();

    const double timestep = 0.0003 * std::sqrt(mass * sigma * sigma / epsilon);
    const double total_time = 100 * std::sqrt(mass * sigma * sigma / epsilon);
    const int num_it = total_time / timestep;
    double e_pot = 0;
    e_pot = lj_neighbors_list(atoms, neighbor_list, epsilon, sigma, lj_cutoff, lj_cutoff_engergy);
    std::cout << "Simulating for at total time of: "<< total_time << "With at timestep of: " << timestep << "\n";
    std::cout << "Resulting in a total of "<< num_it << " itterations.\n";
    std::ofstream traj("milestones/05/traj.xyz");

    for (int i = 0; i < num_it; ++i) {
        verlet_step1(atoms, timestep);
        e_pot = lj_neighbors_list(atoms, neighbor_list, epsilon, sigma, lj_cutoff, lj_cutoff_engergy);
        verlet_step2(atoms, timestep);
        berendsen_thermostat(atoms,.2,timestep,timestep*1000);


        if(i % 1000 == 0){
            std::cout << "Itteration: " << i << "\n";
            std::cout << "E_total: " << atoms.kinetic_engergy() + e_pot << "\n";
            std::cout << "E_kin: " << atoms.kinetic_engergy() << "\n";
            std::cout << "E_pot: " << e_pot<< "\n";
            std::cout << "Temperatur: " << atoms.temperatur() << "\n";
            write_xyz(traj, atoms);

        }
    }
    traj.close();
}