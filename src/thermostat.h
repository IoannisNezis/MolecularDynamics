//
// Created by ianni on 7/5/23.
//

#ifndef MY_MD_CODE_THERMOSTAT_H
#define MY_MD_CODE_THERMOSTAT_H
#include "atoms.h"
void berendsen_thermostat(Atoms &atoms, double temperature, double timestep, double relaxation_time);
void berendsen_thermostat(Atoms &atoms, int nb_local, double temperature, double timestep, double relaxation_time);

#endif //MY_MD_CODE_THERMOSTAT_H
