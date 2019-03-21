#ifndef RETINA_UTILS_H
#define RETINA_UTILS_H

#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include "eigen_types.h"

//namespace utils {
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> gridMats(int nrows, int ncols);
Eigen::MatrixXb circleMask(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius);
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> rotateGrids(double origin[2], Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double degrees);
Eigen::MatrixXb rectMask(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double orient, double width, double height);
//int test_ops();
//}

#endif //RETINA_UTILS_H
