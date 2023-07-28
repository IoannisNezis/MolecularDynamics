//
// Created by ianni on 7/12/23.
//

#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include <gtest/gtest.h>


void calibrate_system(Atoms &atoms, const NeighborList &neighbor_list, const int target_temperature, const int timestep_fs=1) {

    ducastelle(atoms, neighbor_list);
    for (int i = 0; i < 30000; ++i) {
        if (i%500 == 0){
            std::cout << "Temperature: " << atoms.temperatur() << std::endl;
        }
        verlet_step1(atoms, timestep_fs);
        ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, timestep_fs);
        berendsen_thermostat(atoms, target_temperature, timestep_fs, 500);
    }

    std::cout << "Reached Temperature: " << atoms.temperatur() <<" waiting for 5000 fs"<< std::endl;
    for (int i = 0; i < 5000; ++i) {
        verlet_step1(atoms, timestep_fs);
        ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, timestep_fs);
    }
}
void create_lattice(int nx,int ny,int nz,double lattice_constant){
    std::ofstream file("tests/lattice.xyz");
    file << nx*ny*nz << std::endl << std::endl;
    float dx,dy,dz;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                dx = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                dy = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                dz = static_cast <float> (rand()) / static_cast <float> (RAND_MAX) -0.5;
                file << std::setw(2) << "AU" << " "
                     << std::setw(10) << x*lattice_constant << " " << y*lattice_constant << " " << z*lattice_constant
                     << std::setw(13) << dx << " " << dy << " " << dz
                     << std::endl;
            }
        }
    }
    file.close();
}
TEST(TIMESTEP,GOLD){
    create_lattice(10,10,10,3);

    double timestep_fs = 2; // in fempto seconds
    const int total_simulation_time_fs = 300000;
    const int iterations_per_frame = 1000;
    const int frames = total_simulation_time_fs/iterations_per_frame / timestep_fs;
    const double neighbour_cutoff = 5.5;
    double e_pot = 0;
    const int start_temperature = 1000;


    auto [names, positions, velocities]{read_xyz_with_velocities("tests/lattice.xyz")};
    Masses_t masses {names.size()};
    masses.setOnes();
    masses *= 196.966570; // Mass of Gold
    masses /= 0.009649;   // Adjust to fix the Unit of time
    Atoms atoms = Atoms(names, positions, masses);
    NeighborList neighbor_list(neighbour_cutoff);
    neighbor_list.update(atoms);

    std::ofstream traj("traj.xyz");

    ducastelle(atoms, neighbor_list);
    verlet_step1(atoms, timestep_fs);
    e_pot = ducastelle(atoms, neighbor_list);
    verlet_step2(atoms, timestep_fs);

    std::cout << "Calibrating the system to a start temperature of: " << start_temperature << std::endl;
    calibrate_system(atoms,neighbor_list,start_temperature);
    std::cout << "Done : " << atoms.temperatur() << std::endl;
    // Repeat

    const double E_init = e_pot + atoms.kinetic_engergy();


    for (int i = 0; i < frames*iterations_per_frame; ++i) {

        verlet_step1(atoms, timestep_fs);
        e_pot = ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, timestep_fs);
        if(i % iterations_per_frame == 0){
            std::cout << "Timestep: " << timestep_fs <<" fs" << " -->  E_total: " << atoms.kinetic_engergy() + e_pot << std::endl;
            write_xyz(traj, atoms);
            timestep_fs += 0.1;
        }

        // ramp up timestep

    }
    traj.close();
}