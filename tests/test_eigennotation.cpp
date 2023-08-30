//
// Created by ianni on 7/10/23.
//
#include <gtest/gtest.h>
#include <Eigen/Dense>

using namespace Eigen;
TEST(EIGEN_NOTATION, MATH) {

    Eigen::Array3d domain_length = {3, 1, 2};
    double x = domain_length.row(2).sum();
    std::cout << x << std::endl;
}