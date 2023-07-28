//
// Created by ianni on 7/5/23.
//
#include "thermostat.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep, double relaxation_time){
    double lambda = std::sqrt( 1+(temperature/atoms.temperatur()-1)*(timestep/relaxation_time));
    atoms.velocities *= lambda;
}


void berendsen_thermostat(Atoms &atoms, int nb_local, double temperature, double timestep, double relaxation_time){
    double lambda = std::sqrt( 1+(temperature/atoms.local_temperatur(nb_local)-1)*(timestep/relaxation_time));
    atoms.velocities *= lambda;
}
