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
    int numbers[] = {55, 147, 309, 561, 923, 1415, 2057, 2869, 3871, 5083, 6525, 8217, 10179, 12431, 14993, 17885,
                     21127, 24739, 28741};
    for (int i = 0; i < 19; ++i) {
        const double timestep_fs = 1;

        double EAM_cutoff = 5.0;
        // Get MPI info
        int size;
        int rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        bool is_main_rank = rank == 0;

        int magic_number = numbers[i];
        std::ofstream traj;
        std::ofstream data;
        if (is_main_rank) {
            data = std::ofstream("figures/time/time_per_step_size-" + std::to_string(magic_number) + "_ranks-" +
                                 std::to_string(size) + ".dat");
        }
        // Initializing the System
        if (is_main_rank) {
            std::cout << "Input: cluster_" << magic_number << std::endl;
        }

        auto [names, positions]{read_xyz("clusters/cluster_" + std::to_string(magic_number) + ".xyz")};
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


        domain.enable(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);
        NeighborList neighbor_list(EAM_cutoff);
        neighbor_list.update(atoms);
        // Calc the Forces
        ducastelle(atoms, neighbor_list, EAM_cutoff);
        // Writing Simulations Specs

        auto start = std::chrono::high_resolution_clock::now();
        float itterations = 100;
        const double iterations_per_frame = 10;
        for (int i = 0; i <= itterations; ++i) {
            // Write movement data
            if (i % int(iterations_per_frame) == 0 && i != 0) {
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                if (is_main_rank) {
                    data << duration.count() / iterations_per_frame << std::endl;
                }
                start = stop;
            }
            verlet_step1(atoms, timestep_fs);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            verlet_step2(atoms, timestep_fs);
        }
        data.close();
    }
    MPI_Finalize();
}