#ifndef __VERLET_H
#define __VERLET_H

#include <Eigen/Core>
#include "types.h"
#include "atoms.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double dt);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Forces_t &forces, const Masses &masses, double dt);

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double dt);

void verlet_step2(Velocities_t &velocities, const Forces_t &forces,const Masses &masses, double dt);


void verlet_step1(Atoms &atoms, double dt);
void verlet_step2(Atoms &atoms, double dt);
#endif  // __VERLET_H