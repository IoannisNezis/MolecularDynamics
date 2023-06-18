//
// Created by ianni on 5/10/23.
//

#ifndef MY_MD_CODE_TYPES_H
#define MY_MD_CODE_TYPES_H

#include <Eigen/Dense>

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses = Eigen::ArrayXd;
using Names_t = std::vector<std::string>;

#endif //MY_MD_CODE_TYPES_H
