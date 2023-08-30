//
// Created by ianni on 8/1/23.
//

#ifndef PHYSICS_H
#define PHYSICS_H
#include "Eigen/Dense"
const double boltzmann_constant = 8.617333262 * std::pow(10,-5);
double kinetic_energy_to_temperature(const double e_kin, const int nb_atoms);
double temperatur_to_kinetic_energy(const double temp, const int nb_atoms);

#endif //PHYSICS_H
