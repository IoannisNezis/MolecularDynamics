#include "mpi_support.h"
#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"

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


void melt(const int magic_number,
          const double timestep_fs = 3,
          const int iterations_per_frame = 1000,
          const double start_temperature = 500,
          const double deltaT = 10,
          const double relaxation_time_fs = 10000,
          const double measure_time_fs = 10000
) {

    double EAM_cutoff = 5.0;
    double e_pot = 0;
    double e_init;
    bool pre_heat = true;
    const double pre_heat_duration = 2000;
    const double post_pre_heat_relatxation_time = 10000;
    double phase_timer = 0;
    // Get MPI info
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;
    // Initializing the System
    if (is_main_rank) {
        std::cout << "Input: cluster_" << magic_number << std::endl;
    }

    auto [names, positions]{read_xyz("clusters/cluster_" + std::to_string(magic_number) + ".xyz")};
    positions *= 0.978;
    // adjust positions
    double cluster_radius = positions.row(0).maxCoeff();
    double padding = 6;
    positions += cluster_radius + padding / 2;

    Masses_t masses{names.size()};
    masses.setOnes();
    masses *= 196.966570; // Mass of Gold
    masses /= 0.009649;   // Adjust to fix the Unit of time
    Atoms atoms = Atoms(names, positions, masses);
    // Initializing the Domain
    Domain domain(MPI_COMM_WORLD,
                  {cluster_radius * 2 + padding, cluster_radius * 2 + padding, cluster_radius * 2 + padding},
                  {2, 2, 2},
                  {1,1,1});

    // Initializing IO
    std::ofstream traj;
    std::ofstream data;
    int data_count = 0;
    if (is_main_rank) {
        traj = std::ofstream("movement_plots/traj_" + std::to_string(magic_number) + ".xyz");
        write_xyz(traj, atoms);
        data = std::ofstream("figures/" + std::to_string(magic_number) + ".dat");
        data << "#Time E_pot E_kin Temperature" << std::endl;
    }

    NeighborList neighbor_list(EAM_cutoff);
    neighbor_list.update(atoms);
    // Calc the Forces
    e_init = ducastelle(atoms, neighbor_list, EAM_cutoff);
    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);
    neighbor_list.update(atoms);
    // Init and update the neigbors list
    // Start Simulation
    int i = 0;
    while (true) {

        // Write Data out
        if (i % iterations_per_frame == 0) {

            domain.disable(atoms);
            neighbor_list.update(atoms);
            if (is_main_rank) {
                write_xyz(traj, atoms);
                e_pot = ducastelle(atoms, neighbor_list, EAM_cutoff);
//                std::cout << "Frame: " << i/iterations_per_frame << " kelvin" << std::endl;
//                std::cout << "Time: " << i * timestep_fs << " fs" << std::endl;
//                std::cout << "Potential Energy: "  << e_pot  - e_init << std::endl;
//                std::cout << "Kinetic Energy: "  << atoms.kinetic_engergy() << std::endl;
//                std::cout << "Total Energy: " << e_pot +  atoms.kinetic_engergy() - e_init << " eV" << std::endl;
//                std::cout << "Temperature: "  << atoms.temperatur() << " kelvin" << std::endl;
                data << data_count << " " << e_pot - e_init << " " << atoms.kinetic_engergy() << " "
                     << atoms.temperatur() << std::endl;
                data_count++;

            }
            if(is_main_rank)
                std::cout << atoms.temperatur() << std::endl;
            domain.enable(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
        }
        i++;

        // First Verlet Step
        verlet_step1(atoms, timestep_fs);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        // Caclulate the Forces using Ducastelle
        e_pot = ducastelle(atoms, neighbor_list, EAM_cutoff);
        // Second Verlet Step
        // epot reduce
        verlet_step2(atoms, timestep_fs);

        if (pre_heat) {
            if (i * timestep_fs <= pre_heat_duration) {
                berendsen_thermostat(atoms, domain.nb_local(), start_temperature, timestep_fs, 500);
            } else if (i * timestep_fs >= post_pre_heat_relatxation_time + pre_heat_duration) {
                pre_heat = false;
                if (is_main_rank)
                    std::cout << "Pre heating completet" << std::endl;
            }
        } else if (phase_timer >= relaxation_time_fs) {
            phase_timer = 0;
            atoms.velocities *= std::sqrt(1 + deltaT / atoms.temperatur());
        }
        phase_timer += timestep_fs;
    }
    if (is_main_rank) {
        traj.close();
        data.close();
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    melt(3871, 1, 500, 600, 2, 500, 100);
    //melt(6525, 1, 500, 400, 2.5, 100, 5000);
    //melt(8217, 1, 500, 400, 2.5, 100, 5000);
    //melt(10179, 1, 500, 400, 2.5, 100, 5000);
    //melt(28741, 1, 500, 400, 2.5, 200, 5000);

    MPI_Finalize();
}
