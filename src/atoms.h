//
// Created by ianni on 5/17/23.
//

#ifndef MY_MD_CODE_ATOMS_H
#define MY_MD_CODE_ATOMS_H

#include <Eigen/Dense>
#include "types.h"

class Atoms {
public:
    Names_t names;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const Positions_t &p) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }
    Atoms(const Names_t &n, const Positions_t &p) :
            names{n}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }
    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v) :
            names{n},positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

    Atoms(const int nb_atoms) : positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    [[nodiscard]] double kinetic_engergy() const {
        return 0.5 * (1 * velocities.colwise().squaredNorm()).sum();
    }
};

#endif //MY_MD_CODE_ATOMS_H
