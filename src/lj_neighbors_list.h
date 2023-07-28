//
// Created by ianni on 6/7/23.
//

#ifndef MY_MD_CODE_LJ_NEIGHBORS_LIST_H
#define MY_MD_CODE_LJ_NEIGHBORS_LIST_H
#include "atoms.h"
#include "neighbors.h"
#include <Eigen/Dense>
#include "types.h"
double lj_neighbors_list(Atoms &atoms, const NeighborList &neighbor_list, const double epsilon, const double sigma, const double cutoff_radius, const double cutoff_energy);
#endif //MY_MD_CODE_LJ_NEIGHBORS_LIST_H
