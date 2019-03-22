#ifndef RETINA_UTILS_H
#define RETINA_UTILS_H

#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include "eigen_types.h"

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> gridMats(int nrows, int ncols);
Eigen::MatrixXi circleMask(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius);
std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> rotateGrids(double origin[2], Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double degrees);
Eigen::MatrixXi rectMask(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double orient, double width, double height);
void MatrixXiToCSV(std::string fname, Eigen::MatrixXi mat);
void MatrixXdToCSV(std::string fname, Eigen::MatrixXd mat);
#endif //RETINA_UTILS_H
