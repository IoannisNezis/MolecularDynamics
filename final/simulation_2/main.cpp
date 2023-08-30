#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"
#include "physics.h"

void melt(const int magic_number,
          const double timestep_fs = 5,
          const double EAM_cutoff = 7,
          const bool write_movement_data = true,
          const int iterations_per_frame = 1000,
          const double start_temperature = 600,
          const double end_temperature = 1000,
          const double deltaT = 3,
          const double relaxation_time_fs = 5000
) {

    bool pre_heat = true;
    const double pre_heat_duration = 10000;
    const double post_pre_heat_relatxation_time = 10000;
    double phase_timer = 0;
    // Get MPI info
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;
    // Initializing the System
    if (is_main_rank) {
        std::cout << "Input: cluster_" << magic_number << std::endl;
    }

    auto [names, positions]{read_xyz("final/simulation_2/input/cluster_" + std::to_string(magic_number) + ".xyz")};
    positions *= 0.978;
    // adjust positions
    double cluster_radius = positions.row(0).maxCoeff();
    double padding = 15;
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
                  {1, 1, 1});

    // Initializing measurement
    int heat_input_counter = 0;
    int measure_samples = 0;
    double local_kinetic_energy;
    double global_kinetic_energy;
    double global_kinetic_energy_akku;
    double local_potential_energy;
    double global_potential_energy;
    double global_potential_energy_akku;
    double initial_potential_energy;
    double temperature;
    const double deltaQ = temperatur_to_kinetic_energy(deltaT, atoms.nb_atoms());


    // Initializing IO
    std::ofstream traj;
    std::ofstream data;
    if (is_main_rank) {
        if (write_movement_data)
            traj = std::ofstream("final/simulation_2/output/trajectory/" + std::to_string(magic_number) + ".xyz");
        write_xyz(traj, atoms);
        data = std::ofstream("final/simulation_2/output/data/" + std::to_string(magic_number) + ".dat");
        data << "#E_pot E_kin E_opt Temperature" << std::endl;
    }


    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);
    NeighborList neighbor_list(EAM_cutoff);
    neighbor_list.update(atoms);
    // Calc the Forces
    ducastelle(atoms, neighbor_list, EAM_cutoff);
    // Writing Simulations Specs
    if (is_main_rank) {
        std::cout << "Starting Simulation of Cluster: " << magic_number << std::endl
                  << "Simulation Timestep: " << timestep_fs << " fs" << std::endl
                  << "EAM cutoff: " << EAM_cutoff << " Angström" << std::endl
                  << "pre-heat temperature: " << start_temperature << " Kelvin" << std::endl
                  << "The Simulation will end when a Temperature of: " << end_temperature << " Kelvin is reached"
                  << std::endl
                  << "The Heat-Input is: " << deltaT << " Kelvin" << std::endl
                  << "That corresponds to: " << deltaQ << " eV" << std::endl
                  << "There will be " << (write_movement_data ? "a" : "no") << " movement plot" << std::endl
                  << "There will be " << iterations_per_frame << " iterations per Frame" << std::endl;
    }
    int i = 0;
    while (true) {
        // Measure Energies
        if (!pre_heat && phase_timer >= relaxation_time_fs / 2) {
            local_kinetic_energy = atoms.local_kinetic_engergy(domain.nb_local());
            global_kinetic_energy = MPI::allreduce(local_kinetic_energy, MPI_SUM, MPI_COMM_WORLD);
            global_kinetic_energy_akku += global_kinetic_energy;
            local_potential_energy = ducastelle(atoms, domain.nb_local(), neighbor_list, EAM_cutoff);
            global_potential_energy = MPI::allreduce(local_potential_energy, MPI_SUM, MPI_COMM_WORLD);
            global_potential_energy_akku += global_potential_energy;
            measure_samples++;
        }

        // Write movement data
        if (write_movement_data && !pre_heat && i % iterations_per_frame == 0) {
            domain.disable(atoms);
            neighbor_list.update(atoms);
            global_potential_energy = ducastelle(atoms, neighbor_list, EAM_cutoff);
            if (is_main_rank) {
                write_xyz(traj, atoms);
                std::cout << "Frame: " << i / iterations_per_frame << std::endl;
                std::cout << "Temperature: " << int(atoms.temperatur()) << " Kelvin" << std::endl;
                std::cout << "Total Energy: " << atoms.kinetic_engergy() + global_potential_energy << " eV"
                          << std::endl;
            }
            domain.enable(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
        }

        verlet_step1(atoms, timestep_fs);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
        verlet_step2(atoms, timestep_fs);

        if (pre_heat) {
            if (i * timestep_fs <= pre_heat_duration) {
                berendsen_thermostat(atoms, domain.nb_local(), start_temperature, timestep_fs, 1000);
            } else if (i * timestep_fs >= post_pre_heat_relatxation_time + pre_heat_duration) {
                pre_heat = false;
                local_potential_energy = ducastelle(atoms, domain.nb_local(), neighbor_list, EAM_cutoff);
                initial_potential_energy = MPI::allreduce(local_potential_energy, MPI_SUM, MPI_COMM_WORLD);
                if (is_main_rank)
                    std::cout << "Pre heating completet" << std::endl;
            }
        } else {
            if (phase_timer >= relaxation_time_fs) {

                temperature = kinetic_energy_to_temperature(global_kinetic_energy_akku / measure_samples, magic_number);
                if (is_main_rank) {
                    data << global_potential_energy_akku / measure_samples - initial_potential_energy << " "
                         << global_kinetic_energy_akku / measure_samples << " "
                         << heat_input_counter * deltaQ << " "
                         << temperature
                         << std::endl;
                }
                if (end_temperature < temperature) {
                    std::cout << "Goal Temperature reached" << std::endl;
                    break;
                }
                measure_samples = 0;
                global_kinetic_energy_akku = 0;
                global_potential_energy_akku = 0;

                heat_input_counter++;
                phase_timer = 0;
                //atoms.velocities *= std::sqrt(1 + deltaT / atoms.local_temperatur(domain.nb_local()));
            }

            phase_timer += timestep_fs;
        }
        i++;
    }
    if (is_main_rank) {
        traj.close();
        data.close();
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    int numbers[] = {55, 147, 309, 561, 923, 1415, 2057, 2869, 3871, 5083, 6525, 8217, 10179, 12431, 14993, 17885,
                     21127, 24739, 28741};
    //melt(numbers[10],12,true,2000,600,1100,2,2000);
    //melt(numbers[10],5,true,200,500,900,0.5,8000);
    if (argc != 10) {
        if (is_main_rank)
            std::cout << "Invalid arguments amount (" << argc << "), provide cube size (int) as argument." << std::endl;
        MPI_Finalize();
        return 0;
    }
    int magic_number = std::stoi(argv[1]);
    double timestep_fs = std::stoi(argv[2]);
    double EAM_cutoff = std::stoi(argv[3]);
    bool movement_plot = argv[4];
    int iterations_per_frame = std::stoi(argv[5]);
    double start_temperature = std::stoi(argv[6]);
    double end_temperature = std::stoi(argv[7]);
    double deltaT = std::stoi(argv[8]);
    double relaxation_time = std::stoi(argv[9]);
    melt(magic_number, timestep_fs, EAM_cutoff, movement_plot, iterations_per_frame, start_temperature, end_temperature,
         deltaT, relaxation_time);

    MPI_Finalize();
}
