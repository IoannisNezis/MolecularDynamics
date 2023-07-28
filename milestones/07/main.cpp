//
// Created by ianni on 7/5/23.
//
#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include <sstream>

void calibrate_system(Atoms &atoms,
                      NeighborList &neighbor_list,
                      const double target_temperature,
                      const double timestep_fs = 1) {

    ducastelle(atoms, neighbor_list);
    verlet_step1(atoms, timestep_fs);
    neighbor_list.update(atoms);
    ducastelle(atoms, neighbor_list);
    verlet_step2(atoms, timestep_fs);
    for (int i = 0; i < 20000; ++i) {
        if (i % 5000 == 0) {
            std::cout << "Reached Temperature: " << atoms.temperatur() << std::endl;
        }
        verlet_step1(atoms, timestep_fs);
        ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, timestep_fs);
        berendsen_thermostat(atoms, target_temperature, timestep_fs, 3000);
    }

    std::cout << "Reached Temperature: " << atoms.temperatur() << " waiting for 10000 fs" << std::endl;
    for (int i = 0; i < 10000; ++i) {
        verlet_step1(atoms, timestep_fs);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list);
        verlet_step2(atoms, timestep_fs);
    }
}


void melt(  const int magic_number,
              const double timestep_fs = 3,
              const double total_simulation_time_fs = 5000000,
              const int iterations_per_frame = 1000,
              const double neighbour_cutoff = 5.5,
              const double start_temperature = 500,
              const double deltaT = 10,
              const double relaxation_time_fs = 10000,
              const double measure_time_fs = 10000
) {
    double e_pot = 0;
    // heat utility
    bool measure = true; // true means masure, false means realax after heat input
    double measure_buffer = 0;
    double min_temp, max_temp;
    int measure_counter = 0;
    double phase_timer = 0;
    float measure_buffer2 = 0;

    // Initializing the System
    //create_lattice(10  ,10,10,2.8* sigma);
    std::cout << "Input: cluster_" <<magic_number << std::endl;


    auto [names, positions]{read_xyz("clusters/cluster_" + std::to_string(magic_number) + ".xyz")};
    Masses_t masses{names.size()};
    masses.setOnes();
    masses *= 196.966570; // Mass of Gold
    masses /= 0.009649;   // Adjust to fix the Unit of time
    Atoms atoms = Atoms(names, positions, masses);
    NeighborList neighbor_list(neighbour_cutoff);
    neighbor_list.update(atoms);


    std::cout << "Simulating: " << total_simulation_time_fs << "fs (femto seconds)" << std::endl;
    std::cout << "Time step: " << timestep_fs << "fs" << std::endl;
    std::cout << "Simulating " << total_simulation_time_fs / iterations_per_frame / timestep_fs
              << " frames, each frame takes: " << iterations_per_frame << " itterations" << std::endl;
    std::cout << "Resulting in a total of " << total_simulation_time_fs / timestep_fs << " itterations." << std::endl;

    // Initializing IO
    std::ofstream traj("movement_plots/traj_" + std::to_string(magic_number) +".xyz");
    // Temperatur vs Energy
    std::ofstream fig_1("figures/cluster_" + std::to_string(magic_number) + "_energy_temperature.dat");
    fig_1 << "#Total_Energy AVG_Temperature Min_Temperature Max_Temperature" << std::endl;
    // Kinetic energy vs Time
    std::ofstream fig_2("figures/cluster_" + std::to_string(magic_number) + "_time_Ekin.dat");
    fig_2 << "#Time  Kinetic_Energy" << std::endl;

    std::cout << "Calibrating the system to a start temperature of: " << start_temperature << std::endl;
    calibrate_system(atoms, neighbor_list, start_temperature, timestep_fs);
    std::cout << "Done, starting temperature: " << atoms.temperatur() << std::endl;

    e_pot = ducastelle(atoms, neighbor_list);
    const double E_init = e_pot + atoms.kinetic_engergy();
    max_temp = atoms.temperatur();
    min_temp = atoms.temperatur();
    // Repeat
    for (int i = 0; i * timestep_fs < total_simulation_time_fs; ++i) {
        // First Verlet Step
        verlet_step1(atoms, timestep_fs);
        // Caclulate the Forces using Ducastelle
        neighbor_list.update(atoms);
        e_pot = ducastelle(atoms, neighbor_list);
        // Second Verlet Step
        verlet_step2(atoms, timestep_fs);

        if(i % 100 == 0 && i != 0){
            fig_2 << i*timestep_fs << " " << measure_buffer2 / 100 << std::endl;
            measure_buffer2 = 0;
        }else{
            measure_buffer2 += atoms.temperatur();
        }

        phase_timer += timestep_fs;
        if (measure) {
            measure_buffer = measure_buffer + atoms.temperatur();
            measure_counter++;

            min_temp = std::min(min_temp, atoms.temperatur());
            max_temp = std::max(max_temp, atoms.temperatur());

            if (phase_timer >= measure_time_fs) {
                measure_buffer /= measure_counter;
                std::cout << "Measurement completed after " << measure_counter << " samples" << std::endl;
                std::cout << "AVG temperature: " << measure_buffer << " Kelvin" << std::endl;
                std::cout << "Total Energy: " << e_pot + atoms.kinetic_engergy() - E_init << " eV" << std::endl;
                fig_1 <<  e_pot + atoms.kinetic_engergy() - E_init << " " << measure_buffer << " " << min_temp << " "
                      << max_temp << std::endl;
                measure_counter = 0;
                measure_buffer = 0;
                measure = false;
                phase_timer = 0;

                atoms.velocities *=  std::sqrt(1 + deltaT / atoms.temperatur());

                std::cout << "==========================================" << std::endl;
                std::cout << "Frame " << i / iterations_per_frame << std::endl;
                std::cout << "Introducing " << deltaT << " kelvin, waiting " << relaxation_time_fs << " fs" << std::endl;
            }
        } else if (phase_timer > relaxation_time_fs) {
            measure = true;
            phase_timer = 0;
            max_temp = atoms.temperatur();
            min_temp = atoms.temperatur();
            std::cout << "Measuring average temperature for " << measure_time_fs << " fs" << std::endl;
        }

        if (i % iterations_per_frame == 0) {
//            std::cout << "------------------------------------------"<< std::endl;
//            std::cout << "Frame: " << i/iterations_per_frame <<", Itteration: " << i << "\n";
//            std::cout << "E_total: " << atoms.kinetic_engergy() + e_pot << "\n";
//            std::cout << "E_kin: " << atoms.kinetic_engergy() << "\n";
//            std::cout << "E_pot: " << e_pot<< "\n";
//            std::cout << "Temperatur: " << atoms.temperatur() << "\n";
            write_xyz(traj, atoms);
        }
    }
    std::cout << "==========================================" << std::endl;
    traj.close();
    fig_1.close();
    fig_2.close();
}

int main(int argc, char *argv[]) {
    // Configuring the Simulation
//    const double timestep_fs = 3,
//    const int total_simulation_time_fs = 5000000,
//    const int iterations_per_frame = 1000,
//    const double neighbour_cutoff = 5.5,
//    const int start_temperature = 500,
//    const double deltaQ = 0.4,
//    const double relaxation_time_fs = 10000,
//    const int measure_time_fs = 10000


    // melt(55, 1, 1200000,1000,5.5,200,100,10000,5000);
    //melt(147, 1, 2000000,1000,5.5,600,0.5,50000,5000);
    //melt(309, 1, 2000000,1000,5.5,350,0.3,15000,5000);
  //  melt(561, 1, 2500000,1000,5.5,400,0.7,20000,5000);
    melt(561, 1, 2500000,1000,5.5,400,0.7,20000,5000);
    //melt(2057,1,1000000, 500, 5.5,620,20,15000,3000);
}
