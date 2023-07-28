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
    Masses_t masses;

    Atoms(const Positions_t &p) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
    Atoms(const Names_t &n, const Positions_t &p) :
            names{n}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }
    Atoms(const Names_t &n, const Positions_t &p, const Masses_t &m) :
            names{n}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{m} {
        velocities.setZero();
        forces.setZero();
    }
    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v) :
            names{n},positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()}  {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
    }

    Atoms(const int nb_atoms) : positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
    }

    size_t nb_atoms() const {
        return positions.cols();
    }

    double kinetic_engergy() const {
        return 0.5 * ((velocities.colwise().squaredNorm()).rowwise() * masses.transpose()).sum();
    }
    double temperatur() const {
        double boltzmann_constant = 8.617333262 * std::pow(10,-5);
        return 2*kinetic_engergy()/(boltzmann_constant * nb_atoms()*3.0);
    }
    double local_kinetic_engergy(int nb_local) const {
        return 0.5 * ((velocities.colwise().squaredNorm()).rowwise() * masses.transpose()).leftCols(nb_local).sum();
    }
    double local_temperatur(int nb_local) const {
        double boltzmann_constant = 8.617333262 * std::pow(10,-5);
        return 2*local_kinetic_engergy(nb_local)/(boltzmann_constant * nb_local*3.0);
    }
    void resize(const int new_size){
        positions.conservativeResize(3,new_size);
        velocities.conservativeResize(3,new_size);
        forces.conservativeResize(3,new_size);
        masses.conservativeResize(new_size);
        names.resize(new_size);
    }
};

#endif //MY_MD_CODE_ATOMS_H
