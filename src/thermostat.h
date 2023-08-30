//
// Created by ianni on 7/5/23.
//

#ifndef MY_MD_CODE_THERMOSTAT_H
#define MY_MD_CODE_THERMOSTAT_H
#include "atoms.h"

/**
 * Changes the temperature of a system
 *
 * This slowly changes the temperature of a system by rescaling the velocities of each atom.
 * It does that over a timespan tau_relax.
 *
 * @param atoms MD system
 * @param temperature target temperature
 * @param timestep used timestep
 * @param relaxation_time tau_relax
 */
void berendsen_thermostat(Atoms &atoms, double temperature, double timestep, double relaxation_time);

/**
 * Changes the temperature of a system in a decomposed state
 *
 * This is used when using domain decomposition.
 * It uses the local temperature of the system to change the temperature of each subdomain.
 * Note that that is not the same as using the thermostat in a composed state.
 * This function introduces a drift into the system.
 *
 * @param atoms MD system
 * @param temperature target temperature
 * @param timestep used timestep
 * @param relaxation_time tau_relax
 */
void berendsen_thermostat(Atoms &atoms, int nb_local, double temperature, double timestep, double relaxation_time);

#endif //MY_MD_CODE_THERMOSTAT_H
