#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"
#include "physics.h"
#include "lattice.h"
#include <chrono>
#include <iomanip>
#include <filesystem>

int meaure_energy_conservation(int cube_size){
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    const double lattice_constant = 2.68;
    const int padding = 5;
    if (is_main_rank) {
        if (!std::filesystem::exists("final/simulation_1/input/cube_" + std::to_string(cube_size) + ".xyz")) {
            std::cout << "Input does not exist yet, generating new input file..." << std::endl;
            create_lattice(cube_size, cube_size, cube_size, lattice_constant, padding);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto [names, positions, velocitys]{
            read_xyz_with_velocities("final/simulation_1/input/cube_" + std::to_string(cube_size) + ".xyz")};

    Atoms atoms = Atoms(names, positions, velocitys);
    atoms.set_masses(196.966570 / 0.009649);


    const double timestep_fs = 0.5;
    const int iterations_per_frame = 50;
    double inital_temperature = 700;
    double EAM_cutoff = 5.0;


    std::ofstream traj;
    std::ofstream data;
    if (is_main_rank) {
        traj = std::ofstream(
                "final/simulation_1/output/trajectory/cube_" + std::to_string(cube_size) + ".xyz");
        data = std::ofstream(
                "final/output/trajectory/cube_" + std::to_string(cube_size) + ".dat");
        data << "#Strain% Stress Temperature" << std::endl;
    }

    // Initializing the Domain
    Eigen::Array3d
            domain_length = {int(cube_size * lattice_constant) + 2 * padding,
                             int(cube_size * lattice_constant) + 2 * padding,
                             int(cube_size * lattice_constant) + 2 * padding};
    Domain domain(MPI_COMM_WORLD,
                  domain_length,
                  {1, 1, size},
                  {0, 0, 0});

    domain.enable(atoms);
    domain.update_ghosts(atoms, EAM_cutoff * 2);

    NeighborList neighbor_list(EAM_cutoff);

    neighbor_list.update(atoms);

    // Calc inital the Forces
    ducastelle(atoms, neighbor_list, EAM_cutoff);

    for (int i = 0; i < 100000; ++i) {
        if (i % int(iterations_per_frame) == 0) {
            domain.disable(atoms);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            if (is_main_rank) {
                std::cout << atoms.temperatur() << std::endl;
                write_xyz(traj, atoms);

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
    }
    traj.close();
    data.close();
    MPI_Finalize();

}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    if (argc != 2) {
        if (is_main_rank)
            std::cout << "Invalid arguments, provide cube size (int) as argument." << std::endl;
        MPI_Finalize();
        return 0;
    }
    int cube_size = std::stoi(argv[1]);

    MPI_Finalize();
    return 0;
}