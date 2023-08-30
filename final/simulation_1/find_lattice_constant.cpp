#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "domain.h"
#include "physics.h"
#include "lattice.h"

void find_lattice_constant() {
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;
    const int padding = 5;
    const int cube_size = 15;
    std::ofstream traj;
    std::ofstream data;
    if (is_main_rank) {
        traj = std::ofstream(
                "final/simulation_1/output/trajectory/cube_" + std::to_string(cube_size) + ".xyz");
// data = std::ofstream("final/simulation_1/output/data/LatticeConstant_Temperature.dat");
  //      data << "#Lattice_constant Temperature" << std::endl;
    }

    // Loop over lattice constant
    for (double lattice_constant = 2.5; lattice_constant <= 3; lattice_constant += 0.01) {
        if (is_main_rank) {
            std::cout << "Testing lattice constant: " << lattice_constant << std::endl;
            create_lattice(cube_size, cube_size, cube_size, lattice_constant, padding);
        }
        // Make sure that the file is written
        MPI_Barrier(MPI_COMM_WORLD);
        auto [names, positions, velocitys]{
                read_xyz_with_velocities("final/simulation_1/input/cube_" + std::to_string(cube_size) + ".xyz")};
        // Initialize MD system
        Atoms atoms = Atoms(names, positions, velocitys);
        // adjust masses for unit system to make sense
        atoms.set_masses(196.966570 / 0.009649);

        // set timestep & EAM-cutoff
        const double timestep_fs = 8;
        const int iterations_per_frame = 10;
        double EAM_cutoff = 5.0;

        // Initializing the Domain
        // None periodic
        // Split along the z axis
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

        // Initialize Neighbor-list
        NeighborList neighbor_list(EAM_cutoff);
        neighbor_list.update(atoms);

        // Calc initial the Forces
        ducastelle(atoms, neighbor_list, EAM_cutoff);

        // Simulate 8000 fs
        for (int i = 0; i * timestep_fs < 8000; ++i) {

            // Write  trajectory file
            if (i % int(iterations_per_frame) == 0) {
                domain.disable(atoms);
                neighbor_list.update(atoms);
                ducastelle(atoms, neighbor_list, EAM_cutoff);
                if (is_main_rank) {
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
        domain.disable(atoms);
        if (is_main_rank) {
            data << lattice_constant << " " << atoms.temperatur() << std::endl;
            std::cout << "Temperature after 8000 fs: " << atoms.temperatur() << " Kelvin" << std::endl;
        }
    }
    traj.close();
    data.close();
}
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    find_lattice_constant();

    MPI_Finalize();
    return 0;
}