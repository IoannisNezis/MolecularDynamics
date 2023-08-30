#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"
#include "physics.h"
#include <chrono>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    const double strain_per_fs = 0.000005; // Strain rate
    const double max_strain = 40;
    const double timestep_fs = 10;
    const int iterations_per_frame = 500;
    //for (int radius = 10; radius <=20; ++radius) {
    int radius = 10;
    int length = 100;
    double inital_temperature = 700;
    double EAM_cutoff = 5.0;
    // Get MPI info
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    std::ofstream traj;
    std::ofstream data;
    if (is_main_rank) {
        traj = std::ofstream(
                "movement_plots/Radius-" + std::to_string(radius) + "_StrainRate-1_T-" +
                std::to_string(int(inital_temperature)) + ".xyz");
        data = std::ofstream(
                "figures/stress/strain_vs_stress_Radius-" + std::to_string(radius) + "_StrainRate-1_T-" +
                std::to_string(int(inital_temperature)) + ".dat");
        data << "#Strain% Stress Temperature" << std::endl;
        data << "Radius:" << radius << ", "
             << "Length:" << length << ", "
             << "StrainRate:" << strain_per_fs << ", "
             << "InitalTemp:" << inital_temperature << ", "
             << "Timestep:" << timestep_fs << ", "
             << "EAM-cutoff:" << EAM_cutoff << std::endl;
    }
    auto [names, positions]{read_xyz("whiskers/whisker_R" + std::to_string(radius) + ".xyz")};

    Atoms atoms = Atoms(names, positions);
    atoms.set_masses(196.966570 / 0.009649);

    double max_x = positions.row(0).maxCoeff();
    double min_x = positions.row(0).minCoeff();
    double max_y = positions.row(1).maxCoeff();
    double min_y = positions.row(1).minCoeff();
    double max_z = positions.row(2).maxCoeff();

    if (is_main_rank)
        std::cout << "Object: Whisker, Radius: " << radius << ", Length " << length << std::endl
                  << "Strai-Rate: " << strain_per_fs << " strain per fs" << std::endl
                  << "Timestep: " << timestep_fs << " fs" << std::endl
                  << "EAM Cutoff: " << EAM_cutoff << std::endl;


    // Initializing the Domain
    Eigen::Array3d
            domain_length = {int(min_x + max_x), int(min_y + max_y), int(max_z) + 1};
    Domain domain(MPI_COMM_WORLD,
                  domain_length,
                  {1, 1, 8},
                  {0, 0, 1});

    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);

    NeighborList neighbor_list(EAM_cutoff);

    neighbor_list.update(atoms);

    // Calc the Forces
    ducastelle(atoms, neighbor_list, EAM_cutoff);
    // Writing Simulations Specs

    auto start = std::chrono::high_resolution_clock::now();
    float itterations = max_strain / (strain_per_fs * timestep_fs);


    double global_stress, local_stress;
    double global_kineticE, local_kineticE;
    for (int i = 0; i <= 3000; ++i) {
        verlet_step1(atoms, timestep_fs);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
        verlet_step2(atoms, timestep_fs);
        berendsen_thermostat(atoms, domain.nb_local(), inital_temperature, timestep_fs, 2000);
    }
    domain.disable(atoms);
    if (is_main_rank) {
        std::cout << "Temperature: " << atoms.temperatur() << std::endl;
    }
    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);
    neighbor_list.update(atoms);
    for (int i = 0; i <= itterations; ++i) {
        // Write movement data
        if (i % int(iterations_per_frame) == 0 && i != 0) {
            domain.disable(atoms);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            if (is_main_rank) {
                write_xyz(traj, atoms);
                //std::cout << global_stress << std::endl;
                data << 100 * (i * (strain_per_fs * timestep_fs)) / max_z << " " << global_stress / iterations_per_frame
                     << " " << kinetic_energy_to_temperature(global_kineticE/iterations_per_frame, atoms.nb_atoms())
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
        verlet_step1(atoms, timestep_fs);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
        verlet_step2(atoms, timestep_fs);

        //Measure Stress
        local_stress =
                1 / (2 * domain_length.prod()) * atoms.partial_stress.leftCols(domain.nb_local()).rowwise().sum().sum();
        global_stress += MPI::allreduce(local_stress, MPI_SUM, MPI_COMM_WORLD);
        // Measure Temperature
        local_kineticE = atoms.local_kinetic_engergy(domain.nb_local());
        global_kineticE += MPI::allreduce(local_kineticE, MPI_SUM, MPI_COMM_WORLD);
        //strech the cluster
        domain_length.row(2) += strain_per_fs * timestep_fs;
        domain.scale(atoms, domain_length);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        neighbor_list.update(atoms);
        ducastelle(atoms, neighbor_list, EAM_cutoff);
    }
    //   }
    MPI_Finalize();
}