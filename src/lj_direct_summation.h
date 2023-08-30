//
// Created by ianni on 5/17/23.
//

#ifndef MY_MD_CODE_LS_DIRECT_SUMMATION_H
#define MY_MD_CODE_LS_DIRECT_SUMMATION_H
#include "atoms.h"

/**
 * Calculates the Potential energy for a system and update the forces.
 *
 * Calculates the potential energy using the Lenard jones potential.
 * The forces of each atom is calulated using the derivation.
 *
 * @param atoms MD system
 * @param epsilon epsilon: how deep is the valley
 * @param sigma sigma: When does the repulsion start -> radius of atoms
 * @return
 */
double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

double lj_potential(double distance, double epsilon, double sigma);
#endif //MY_MD_CODE_LS_DIRECT_SUMMATION_H
