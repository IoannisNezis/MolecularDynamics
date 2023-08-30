#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Core>
#include "types.h"
#include "atoms.h"

/**
 * First verlet step
 *
 * This claculates v(t+ Delta t/2) and updates the position.
 *
 * @param atoms MD system
 * @param dt timestep Delta t
 */
void verlet_step1(Atoms &atoms, double dt);

/**
 * Second verlet step
 *
 * This calculates v(t+Delta t).
 * Its requred to update the forces between verlet step 1 and 2.
 *
 * @param atoms MD system
 * @param dt timestep Delta t
 */
void verlet_step2(Atoms &atoms, double dt);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, const Masses_t &masses, double dt);

void verlet_step2(Velocities_t &velocities, const Forces_t &forces,const Masses_t &masses, double dt);


void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double dt);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double dt);

#endif  // __VERLET_H