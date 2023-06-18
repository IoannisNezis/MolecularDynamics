//
// Created by ianni on 5/17/23.
//
#include <iostream>
#include "lj_direct_summation.h"
#include <Eigen/Dense>

double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    atoms.forces.setZero();
    auto N = atoms.nb_atoms();
    double Epot = 0;
    Eigen::Vector3d distance;
    double forceScale;
    double temp;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            distance = atoms.positions.col(j) - atoms.positions.col(i);
            temp = std::pow(sigma / distance.norm(), 6);
            Epot += 4 * epsilon * (std::pow(temp, 2) - temp);
            forceScale = -24 * epsilon / distance.norm() * (2 * std::pow(temp, 2) - temp);
            atoms.forces.col(i)   += forceScale * distance.normalized().array();
            atoms.forces.col(j) -= forceScale * distance.normalized().array();
        }
    }
    return Epot;
}


