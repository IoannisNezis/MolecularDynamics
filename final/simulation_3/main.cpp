#include <iostream>
#include <iomanip>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"
#include "physics.h"
#include <chrono>
#include <filesystem>


void pull(double timestep_fs, double EAM_cutoff, int radius, double strain_rate, double max_strain,
          double inital_temperature) {
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    int length = 100;
    // Calculate strech per fs based on strain rate
    double strech_per_fs = length * strain_rate / 100.0;
    int iterations_per_frame = 100;
    float iterations = max_strain / (strain_rate * timestep_fs);
    std::ofstream traj;
    std::ofstream data;
    if (is_main_rank) {
        traj = std::ofstream(
                "final/simulation_3/output/trajectory/Radius-" + std::to_string(radius) + "_StrainRate-" +
                std::to_string(strain_rate) + "_T-" +
                std::to_string(int(inital_temperature)) + ".xyz");
        data = std::ofstream(
                "final/simulation_3/output/data/strain_vs_stress_Radius-" + std::to_string(radius) + "_StrainRate-" +
                std::to_string(strain_rate) + "_T-" +
                std::to_string(int(inital_temperature)) + ".dat");
        data << "#Strain% Stress Temperature" << std::endl;
        data << "Radius:" << radius << ", "
             << "Length:" << length << ", "
             << "StrainRate:" << strain_rate << ", "
             << "InitalTemp:" << inital_temperature << ", "
             << "Timestep:" << timestep_fs << ", "
             << "EAM-cutoff:" << EAM_cutoff << std::endl;
    }
    // Initialize MD system
    auto [names, positions]{read_xyz("final/simulation_3/input/whisker_R" + std::to_string(radius) + ".xyz")};
    Atoms atoms = Atoms(names, positions);
    // If there already exists a inputfile use it
    if (std::filesystem::exists("final/simulation_3/input/Radius-" + std::to_string(radius) + "_T-" +
                                std::to_string(int(inital_temperature)) + ".xyz")) {
        auto [names, positions, velocities]{read_xyz_with_velocities(
                "final/simulation_3/input/Radius-" + std::to_string(radius) + "_T-" +
                std::to_string(int(inital_temperature)) + ".xyz")};
        atoms = Atoms(names, positions, velocities);
    }
    // Adjust Masses
    atoms.set_masses(196.966570 / 0.009649);

    double max_x = positions.row(0).maxCoeff();
    double min_x = positions.row(0).minCoeff();
    double max_y = positions.row(1).maxCoeff();
    double min_y = positions.row(1).minCoeff();
    double max_z = positions.row(2).maxCoeff();

    if (is_main_rank)
        std::cout << "Object: Whisker, Radius: " << radius << ", Length " << length << std::endl
                  << "Strai-Rate: " << strain_rate << " strain per fs" << std::endl
                  << "Timestep: " << timestep_fs << " fs" << std::endl
                  << "Iterations: " << iterations << std::endl
                  << "EAM Cutoff: " << EAM_cutoff << std::endl;


    // Initializing the Domain
    // Periodic in z direction
    // split subdomains in z direction
    Eigen::Array3d
            domain_length = {int(min_x + max_x), int(min_y + max_y), int(max_z) + 1};
    Domain domain(MPI_COMM_WORLD,
                  domain_length,
                  {1, 1, size},
                  {0, 0, 1});

    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);

    NeighborList neighbor_list(EAM_cutoff);

    neighbor_list.update(atoms);

    // Calc the Forces
    ducastelle(atoms, neighbor_list, EAM_cutoff);
    // Writing Simulations Specs

    auto start = std::chrono::high_resolution_clock::now();


    double global_stress, local_stress;
    double global_kineticE, local_kineticE;

    // If such a input file does not exist yet, heat up and save
    if (!std::filesystem::exists("final/simulation_3/input/Radius-" + std::to_string(radius) + "_T-" +
                                 std::to_string(int(inital_temperature)) + ".xyz")) {
        if(is_main_rank)
            std::cout << "Pre heating to: " << inital_temperature <<" Kelvin" << std::endl;
        for (int i = 0; i <= 3000; ++i) {
            verlet_step1(atoms, timestep_fs);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            verlet_step2(atoms, timestep_fs);
            berendsen_thermostat(atoms, domain.nb_local(), inital_temperature, timestep_fs, 2000);
        }
        for (int i = 0; i <= 5000; ++i) {
            verlet_step1(atoms, timestep_fs);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            verlet_step2(atoms, timestep_fs);
        }
        domain.disable(atoms);
        if (is_main_rank) {
            std::ofstream pre_heat = std::ofstream(
                    "final/simulation_3/input/Radius-" + std::to_string(radius) + "_T-" +
                    std::to_string(int(inital_temperature)) + ".xyz");
            write_xyz(pre_heat, atoms);
            pre_heat.close();
            std::cout << "Temperature: " << atoms.temperatur() << std::endl;
        }
        domain.enable(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
    }

    // Core simulation loop
    for (int i = 0; i <= iterations; ++i) {
        // Write data
        if (i % int(iterations_per_frame) == 0) {
            domain.disable(atoms);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            if (is_main_rank) {
                write_xyz(traj, atoms);
                std::cout << '\r' << "Progress: "
                          << std::setw(4) << 100 * i * (strain_rate * timestep_fs) / max_strain
                          << std::setw(1) << "%     " << std::flush;
                //std::cout << "Strain: " << i * strain_rate * timestep_fs << "%/" << max_strain << "%" << std::endl;
                //std::cout << "Stress: " << global_stress / iterations_per_frame << std::endl;
                data << i * (strain_rate * timestep_fs) << " " << global_stress / iterations_per_frame
                     << " " << kinetic_energy_to_temperature(global_kineticE / iterations_per_frame, atoms.nb_atoms())
                     << std::endl;
                global_stress = 0;
                global_kineticE = 0;
            }
            start = stop;
            domain.enable(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
        }
        // Simulation step
        verlet_step1(atoms, timestep_fs);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
        verlet_step2(atoms, timestep_fs);

        //Measure Stress
        local_stress = 1 / (domain_length.row(2).sum()) *
                       (atoms.partial_stress.leftCols(domain.nb_local()).rowwise().sum().sum());
        global_stress += MPI::allreduce(local_stress, MPI_SUM, MPI_COMM_WORLD);
        // Measure Temperature
        local_kineticE = atoms.local_kinetic_engergy(domain.nb_local());
        global_kineticE += MPI::allreduce(local_kineticE, MPI_SUM, MPI_COMM_WORLD);
        //strech the cluster
        domain_length.row(2) += strech_per_fs * timestep_fs;
        domain.scale(atoms, domain_length);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    double timestep_fs = std::atof(argv[1]);
    double EAM_cutoff = std::atof(argv[2]);
    int radius = std::stoi(argv[3]);
    double strain_rate = std::atof(argv[4]);
    double max_strain = std::atof(argv[5]);
    double initial_temperature = std::atof(argv[6]);

    pull(timestep_fs, EAM_cutoff, radius, strain_rate, max_strain, initial_temperature);
    // Get MPI info

    //   }
    MPI_Finalize();
}