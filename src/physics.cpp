//
// Created by ianni on 8/1/23.
//
#include "physics.h"

double kinetic_energy_to_temperature(const double e_kin, const int nb_atoms){
    return 2*e_kin/(boltzmann_constant * nb_atoms*3.0);
}
double temperatur_to_kinetic_energy(const double temp, const int nb_atoms){
    return temp*boltzmann_constant*nb_atoms*(3/2.0);
}
