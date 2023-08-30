//
// Created by ianni on 8/22/23.
//


#ifndef MY_MD_CODE_LATTICE_H
#define MY_MD_CODE_LATTICE_H
#include <fstream>
#include <iomanip>
/**
 * Creates a cube of golf atoms
 *
 * @param nx number of atoms in the x direction
 * @param ny number of atoms in the y direction
 * @param nz number of atoms in the z direction
 * @param lattice_constant distance between atoms
 * @param padding
 */
void create_lattice(int nx, int ny, int nz, double lattice_constant, int padding) {
    std::ofstream file("final/simulation_1/input/cube_" + std::to_string(nx) + ".xyz");
    file << nx * ny * nz << std::endl << std::endl;
    float dx, dy, dz;
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                // Create random velocities
                dx = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5) * 0.0001;
                dy = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5) * 0.0001;
                dz = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 0.5) * 0.0001;
                file << std::setw(2) << "Au" << " "
                     << std::setw(10) << x * lattice_constant + padding << " " << y * lattice_constant + padding << " "
                     << z * lattice_constant + padding
                     << std::setw(13) << dx << " " << dy << " " << dx
                     << std::endl;
            }
        }
    }
    file.close();
}
#endif //MY_MD_CODE_LATTICE_H
