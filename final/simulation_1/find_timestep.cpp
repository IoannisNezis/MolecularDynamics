#include <iostream>
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include "verlet.h"
#include "thermostat.h"
#include "domain.h"
#include "physics.h"
#include <chrono>
#include <iomanip>
#include <filesystem>
#include "lattice.h"3

void find_timestep() {
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;
    const int padding = 5;
    const int cube_size = 10;
    const double lattice_constant = 2.7;
    const int iterations_per_frame = 100;
    double EAM_cutoff = 7.0;
    std::ofstream traj;
    std::ofstream data1,data2;
    if (is_main_rank) {
        traj = std::ofstream(
                "final/simulation_1/output/trajectory/cube_" + std::to_string(cube_size) + ".xyz");
        data1 = std::ofstream(   "final/simulation_1/output/data/Timestep_Energy.dat");
    }

    for (int timestep_fs = 1; timestep_fs <= 30; timestep_fs += 1) {


        if (is_main_rank) {
            std::cout << "Testing timestep: " << timestep_fs << " fs" << std::endl;
            create_lattice(cube_size, cube_size, cube_size, lattice_constant, padding);

            data2 = std::ofstream("final/simulation_1/output/data/TotalEnerg_"+ std::to_string(timestep_fs)+"_.dat");
            data2 << "#Timestep Energy" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        auto [names, positions, velocitys]{
                read_xyz_with_velocities("final/simulation_1/input/cube_" + std::to_string(cube_size) + ".xyz")};

        Atoms atoms = Atoms(names, positions, velocitys);
        atoms.set_masses(196.966570 / 0.009649);


        double e_pot = 0;
        double e_max = -1000, e_min=1000;
        double local_e_kinetic = 0;
        double global_e_kinetic = 0;
        double local_e_pot = 0;
        double global_e_pot = 0;


        // Initializing the Domain
        Eigen::Array3d
                domain_length = {int(cube_size * lattice_constant) + 2 * padding,
                                 int(cube_size * lattice_constant) + 2 * padding,
                                 int(cube_size * lattice_constant) + 2 * padding};
        Domain domain(MPI_COMM_WORLD,
                      domain_length,
                      {1, 1, size},
                      {1, 1, 1});

        domain.enable(atoms);
        domain.update_ghosts(atoms, EAM_cutoff * 2);

        NeighborList neighbor_list(EAM_cutoff);

        neighbor_list.update(atoms);

        // Calc inital the Forces
        ducastelle(atoms, neighbor_list, EAM_cutoff);

        for (int i = 0; i <100; ++i) {
            verlet_step1(atoms, timestep_fs);
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, EAM_cutoff * 2);
            neighbor_list.update(atoms);
            ducastelle(atoms, neighbor_list, EAM_cutoff);
            verlet_step2(atoms, timestep_fs);

            berendsen_thermostat(atoms, domain.nb_local(), 1000, timestep_fs, timestep_fs*100);
        }
        local_e_kinetic = ducastelle(atoms,domain.nb_local(),neighbor_list,EAM_cutoff);
        const double e_pot_init = MPI::allreduce(local_e_kinetic, MPI_SUM, MPI_COMM_WORLD);

        for (int i = 0; i < 1000; ++i) {
            if (i % int(iterations_per_frame) == 0) {
                domain.disable(atoms);
                neighbor_list.update(atoms);
                e_pot = ducastelle(atoms, neighbor_list, EAM_cutoff);
                if (is_main_rank) {
                    write_xyz(traj, atoms);
                    //std::cout << "Total energy: " << atoms.kinetic_engergy() + e_pot -e_pot_init << std::endl;
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
            local_e_pot = ducastelle(atoms,domain.nb_local(), neighbor_list, EAM_cutoff);
            verlet_step2(atoms, timestep_fs);

            local_e_kinetic = atoms.local_kinetic_engergy(domain.nb_local());
            global_e_kinetic = MPI::allreduce(local_e_kinetic, MPI_SUM, MPI_COMM_WORLD);
            global_e_pot = MPI::allreduce(local_e_pot, MPI_SUM, MPI_COMM_WORLD);
            if(is_main_rank){
                e_max = std::max(e_max, global_e_kinetic + global_e_pot - e_pot_init);
                e_min = std::min(e_min, global_e_kinetic + global_e_pot - e_pot_init);
                data2 << global_e_kinetic + global_e_pot - e_pot_init<< std::endl;
            }
        }
        domain.disable(atoms);
        if (is_main_rank) {
            data1 << timestep_fs << " " << e_max-e_min << std::endl;
            std::cout << "Energy fluctuation: " << e_max-e_min << " eV" << std::endl;
        }
    }
    traj.close();
    //data.close();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bool is_main_rank = rank == 0;

    find_timestep();

    MPI_Finalize();
    return 0;
}