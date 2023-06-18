//
// Created by ianni on 5/3/23.
//
#include "verlet.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  const double fx, const double fy, double fz, double dt) {
    vx += 0.5 * fx * dt; // mass is 1
    vy += 0.5 * fx * dt;
    vz += 0.5 * fx * dt;

    x += vx * dt;
    y += vy * dt;
    z += vz * dt;

}

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, const Masses &masses, double dt) {
    velocities += (0.5 * forces * dt).rowwise() / masses.transpose();
    positions += velocities * dt;
}

void verlet_step1(Atoms &atoms, double dt){
    atoms.velocities += 0.5 * (atoms.forces/1) * dt;
    atoms.positions  += atoms.velocities * dt;
}
void verlet_step2(Atoms &atoms, double dt){
    atoms.velocities += 0.5 * (atoms.forces/1) * dt;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double dt) {
    vx += 0.5 * fx * dt; // mass is 1
    vy += 0.5 * fx * dt;
    vz += 0.5 * fx * dt;
}

void verlet_step2(Velocities_t &velocities, const Forces_t &forces, const Masses &masses, double dt) {
    velocities += (0.5*forces*dt).rowwise()/masses.transpose();
}