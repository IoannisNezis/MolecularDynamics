//
// Created by ianni on 5/17/23.
//

#ifndef MY_MD_CODE_ATOMS_H
#define MY_MD_CODE_ATOMS_H


#include <Eigen/Dense>
#include "types.h"
#include "physics.h"

class Atoms {
public:
    Names_t names;
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses_t masses;
    Eigen::Array3Xd partial_stress;

    Atoms(const Positions_t &p) :
            positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        partial_stress.setZero();
    }

    Atoms(const Names_t &n, const Positions_t &p) :
            names{n}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{p.cols()} {
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        partial_stress.setZero();
    }

    Atoms(const Names_t &n, const Positions_t &p, const Masses_t &m) :
            names{n}, positions{p}, velocities{3, p.cols()}, forces{3, p.cols()}, masses{m} {
        velocities.setZero();
        forces.setZero();
        partial_stress.setZero();
    }

    Atoms(const Names_t &n, const Positions_t &p, const Velocities_t &v) :
            names{n}, positions{p}, velocities{v}, forces{3, p.cols()}, masses{p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
        masses.setOnes();
        partial_stress.setZero();
    }

    Atoms(const int nb_atoms) : positions{3, nb_atoms}, velocities{3, nb_atoms}, forces{3, nb_atoms}, masses{nb_atoms} {
        positions.setZero();
        velocities.setZero();
        forces.setZero();
        masses.setOnes();
        partial_stress.setZero();
    }

    /**
    * Calculates the number of atoms
    *
    * @return number of atoms
     */
    size_t nb_atoms() const {

        return positions.cols();
    }

    /**
     * Sets the mass of each atom
     *
     * @param mass
     */
    void set_masses(const double mass) {

        masses.setOnes();
        masses *= mass;
    }

    /**
     * Calculates the kinetic energy in the system
     *
     * This is done with the formula:
     * E_kin = 1/2 * m * v^2
     *
     * @return the kinetic energy
     */
    double kinetic_engergy() const {
        return 0.5 * ((velocities.colwise().squaredNorm()).rowwise() * masses.transpose()).sum();
    }

    /**
     * Calculates the temperature of the system
     *
     * This is done with the formula:
     * T = 2/3 * E_kin/(k_B * nb_atoms)
     *
     * @return temperature of the system
     */
    double temperatur() const {
        return kinetic_energy_to_temperature(kinetic_engergy(), nb_atoms());
    }

    /**
     * Calculates the kinetic energy in the system for the first nb_local atoms
     *
     * This is usefull when using domain decomposition.
     * This function ignores the ghost atoms
     *
     * @param nb_local number of atoms in the subdomain
     * @return the local kinetic energy
     */
    double local_kinetic_engergy(int nb_local) const {
        return 0.5 * ((velocities.colwise().squaredNorm()).rowwise() * masses.transpose()).leftCols(nb_local).sum();
    }

    /**
     * Calculates the temperature in the system for the first nb_local atoms
     *
     * This is used when using domain decomposition.
     * This function ignores the gost atoms
     *
     * @param nb_local number of atoms in the subdomain
     * @return local temperature
     */
    double local_temperatur(int nb_local) const {
        return kinetic_energy_to_temperature(local_kinetic_engergy(nb_local), nb_local);
    }

    /**
     * Changes the number of atoms in the system.
     *
     * When using domain decomposition, you might want to add or remove atoms from a subdomain.
     * This funktion is used to change the amount of atoms in a subdomain.
     *
     * @param new_size new number of atoms
     */
    void resize(const int new_size) {
        positions.conservativeResize(3, new_size);
        velocities.conservativeResize(3, new_size);
        forces.conservativeResize(3, new_size);
        partial_stress.conservativeResize(3, new_size);
        masses.conservativeResize(new_size);
        names.resize(new_size);
    }
};

#endif //MY_MD_CODE_ATOMS_H
