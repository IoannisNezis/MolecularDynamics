//
// Created by ianni on 5/17/23.
//

#ifndef MY_MD_CODE_LS_DIRECT_SUMMATION_H
#define MY_MD_CODE_LS_DIRECT_SUMMATION_H
#include "atoms.h"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0);

double lj_potential(double distance, double epsilon, double sigma);
#endif //MY_MD_CODE_LS_DIRECT_SUMMATION_H
