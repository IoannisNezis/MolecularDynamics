//
// Created by ianni on 5/3/23.
//
#include <gtest/gtest.h>
#include "verlet.h"
#include "types.h"

void compare(const Eigen::Array3Xd &real, const Eigen::Array3Xd &expected) {
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < real.cols(); ++col) {
            EXPECT_NEAR(real(row, col), expected(row, col), 1e-10);
        }
    }
}

void get_expected_pos_and_vel(Positions_t &exp_pos, Positions_t &init_pos, Velocities_t &exp_vel, Velocities_t &init_vel,
                              const Forces_t &forces, const Masses &masses, const double time) {
    exp_pos = 0.5 * (forces.rowwise() / masses.transpose()) * std::pow(time, 2) + init_vel * time + init_pos;
    exp_vel = (forces.rowwise() / masses.transpose()) * time + init_vel;
}

TEST(VERLETTEST, BasicAssertions) {
    Positions_t a{{1, 2},
                  {1, 2},
                  {1, 2}};

    Masses b{{2},
             {2}};
    int nb_steps = 10;
    int nb_atoms = 2;
    double dt = 1;
    // 1. populate with setRandom() Eigen funktion
    Positions_t positions(3, nb_atoms);
    positions.setRandom();
    Velocities_t velocitys(3, nb_atoms);
    velocitys.setRandom();
    Forces_t forces(3, nb_atoms);
    forces.setRandom();
    Masses masses(nb_atoms);
    masses.setRandom();
    // 2. Initial positions needs to be saved for testting
    Positions_t positions_initial = {positions};
    Positions_t velocitys_initial = {velocitys};
    // 3, Perform Verlet algorithm and Test
    Positions_t expected_pos(3, nb_atoms);
    Positions_t expected_vel(3, nb_atoms);
    for (int i = 0; i < nb_steps; ++i) {

        verlet_step1(positions, velocitys, forces, masses, dt);
        // Compute Forces
        verlet_step2(velocitys, forces, masses, dt);
        get_expected_pos_and_vel(expected_pos, positions_initial, expected_vel, velocitys_initial, forces, masses,
                                 (i+1) * dt);
        compare(positions, expected_pos);
        compare(velocitys, expected_vel);
    }


}
