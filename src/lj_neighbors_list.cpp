//
// Created by ianni on 6/7/23.
//
#include "lj_neighbors_list.h"

double lj_neighbors_list(Atoms &atoms, const NeighborList &neighbor_list, const double epsilon, const double sigma, const double cutoff_radius, const double cutoff_energy) {

    atoms.forces.setZero();
    double Epot = 0;
    Eigen::Vector3d  offset;
    double forceScale;
    double temp;
    for (auto[i, j]: neighbor_list) {
        if (i < j) {
            offset = atoms.positions.col(j) - atoms.positions.col(i);
            if (offset.norm() <= cutoff_radius){
                temp = std::pow(sigma / offset.norm(), 6);
                Epot += 4 * epsilon * (std::pow(temp, 2) - temp) - cutoff_energy;
                forceScale = -24 * epsilon / offset.norm() * (2 * std::pow(temp, 2) - temp);
                atoms.forces.col(i)   += forceScale * offset.normalized().array();
                atoms.forces.col(j) -= forceScale * offset.normalized().array();
            }
        }
    }

    return Epot;
}


