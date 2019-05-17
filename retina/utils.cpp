#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "type_defs.h"
#include "utils.h"

using namespace Eigen;


std::tuple<MatrixXd, MatrixXd> gridMats(int nrows, int ncols) {
    MatrixXd xgrid = MatrixXd::Ones(ncols, nrows); // double
    MatrixXd ygrid = MatrixXd::Ones(ncols, nrows); // double

    // range from 0 to dim_size-1
    VectorXd xgridvec = VectorXd::LinSpaced(ncols, 0, ncols-1);
    VectorXd ygridvec = VectorXd::LinSpaced(nrows, 0, nrows-1);

    // broadcast the range vectors over rows and columns respectively
    xgrid = xgrid.array().rowwise() * xgridvec.transpose().array();
    ygrid = ygrid.array().colwise() * ygridvec.array();

    // package for output
    std::tuple<MatrixXd, MatrixXd> out = std::make_tuple(xgrid, ygrid);
    return out;
}

std::tuple<VectorXd, VectorXd> gridVecs(int nrows, int ncols) {
    // range from 0 to dim_size-1
    VectorXd xvec = VectorXd::LinSpaced(ncols, 0, ncols-1);
    VectorXd yvec = VectorXd::LinSpaced(nrows, 0, nrows-1);

    // package for output
    std::tuple<VectorXd, VectorXd> out = std::make_tuple(xvec, yvec);
    return out;
}

/* Take X and Y grid matrices and use them to calculate the distance from the
 * origin to every element of the matrix. Return a boolean mask identifying
 * all elements that are within a maximum radius (a circle shape).
 */
MatrixXi circleMask(Eigen::VectorXd xgrid, Eigen::VectorXd ygrid, Eigen::VectorXd xOnes,
                    Eigen::VectorXd yOnes, double origin[2], double radius) {
    // squared euclidean distance (not taking sqrt, square the radius instead)
    MatrixXd rgrid = (xgrid.array() - origin[0]).square().matrix()*yOnes.transpose() + xOnes*(ygrid.array() - origin[1]).square().matrix().transpose();
    // convert to boolean based on distance from origin vs radius of desired circle
    MatrixXi mask = (rgrid.array() <= pow(radius, 2)).cast<int>();
    return mask;
}

// Rotate X and Y grid matrices clockwise by given degrees, around an origin.
std::tuple<VectorXd, VectorXd> rotateGrids(double origin[2], VectorXd xgrid, VectorXd ygrid, double degrees) {
    // convert orientation angle to radians
    double theta = deg2rad(degrees);

    // rotate x and y grids around origin
    xgrid = xgrid.array() - origin[0];
    ygrid = ygrid.array() - origin[1];
    MatrixXd x_rot = origin[0] + (cos(theta)*xgrid - sin(theta)*ygrid).array();
    MatrixXd y_rot = origin[1] + (sin(theta)*xgrid + cos(theta)*ygrid).array();

    std::tuple<MatrixXd, MatrixXd> out = std::make_tuple(x_rot, y_rot);
    return out;
}

/* Take X and Y grid matrices and draw a rectangular boolean mask. Rectangle is defined by WxH dims
 * centred on an origin, and oriented in an arbitrary angle.
 */
MatrixXi rectMask(MatrixXd xgrid, MatrixXd ygrid, double origin[2], double orient, double width, double height) {
    MatrixXi mask; // integer

    // rotate coordinates according to orientation of desired rectangle mask
    auto [xgrid_rot, ygrid_rot] = rotateGrids(origin, xgrid, ygrid, orient);
    // convert to boolean based on distances from origin on rotated x and y planes
    mask = (
                ((xgrid_rot.array() - origin[0]).abs() <= width/2).matrix().cast<int>()
                * ((ygrid_rot.array() - origin[1]).abs() <= height/2).matrix().cast<int>().transpose()
            );
    return mask;
}

void MatrixXiToCSV(std::string fname, MatrixXi mat) {
    // CSV format described in eigen_types.h
    std::ofstream file(fname.c_str());
    file << mat.format(CSVFormat);
    file.close();
}

void MatrixXdToCSV(std::string fname, MatrixXd mat) {
    // CSV format described in eigen_types.h
    std::ofstream file(fname.c_str());
    file << mat.format(CSVFormat);
    file.close();
}

double deg2rad(double degrees) {
    return degrees * 3.14159265359/180;
}

