//
// Created by ianni on 8/1/23.
//

#ifndef PHYSICS_H
#define PHYSICS_H
#include "Eigen/Dense"

/**
 * Boltzman constant for our unit system
 */
const double boltzmann_constant = 8.617333262 * std::pow(10,-5);

/**
 * Converts kinetic energy to temperature
 *
 * @param e_kin kinetic energy
 * @param nb_atoms numbe of atoms
 * @return temperature
 */
double kinetic_energy_to_temperature(const double e_kin, const int nb_atoms);

/**
 * Converts temperature to kinetic energy
 *
 * @param temp temperature
 * @param nb_atoms number of atoms
 * @return kinetic energy
 */
double temperatur_to_kinetic_energy(const double temp, const int nb_atoms);

#endif //PHYSICS_H
