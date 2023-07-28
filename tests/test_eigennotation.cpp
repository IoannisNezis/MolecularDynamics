//
// Created by ianni on 7/10/23.
//
#include <gtest/gtest.h>
#include <Eigen/Dense>
using namespace Eigen;
TEST(EIGEN_NOTATION, MATH) {
    Eigen::Array2Xd positions(2, 4);
    positions << 1, 1, 1, 1,
                 1, 2, 3, 4;
    positions += 16;
    std::cout << positions << std::endl;
    auto x = ArrayXd::LinSpaced(10, -1 / 2, 1 / 2);

    ArrayXd u(10);

    // Initial conditions
    u = sin(4 * M_PI * x / 1) * exp(-abs(x));
    std::cout << positions.leftCols(2) << std::endl;
}