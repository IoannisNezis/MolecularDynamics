//
// Created by ianni on 5/17/23.
//
#include <iostream>
#include "lj_direct_summation.h"
#include <Eigen/Dense>

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    // reset the forces
    atoms.forces.setZero();
    auto N = atoms.nb_atoms();
    double Epot = 0;
    Eigen::Vector3d distance;
    double forceScale;
    double temp;
    // for each pair:
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            // compute the distance
            distance = atoms.positions.col(j) - atoms.positions.col(i);
            // temporary term: (sigma/r)^6
            temp = std::pow(sigma / distance.norm(), 6);
            // Lenard jones potential
            Epot += 4 * epsilon * (std::pow(temp, 2) - temp);
            // derivation of the Lenard Jones potential
            forceScale = -24 * epsilon / distance.norm() * (2 * std::pow(temp, 2) - temp);
            // update the force for both atoms
            atoms.forces.col(i)   += forceScale * distance.normalized().array();
            atoms.forces.col(j) -= forceScale * distance.normalized().array();
        }
    }
    return Epot;
}


