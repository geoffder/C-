#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include "eigen_types.h"
#include "utils.h"

using namespace Eigen;


std::tuple<MatrixXd, MatrixXd> gridMats(int nrows, int ncols){
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


/* Take X and Y grid matrices and use them to calculate the distance from the
 * origin to every element of the matrix. Return a boolean mask identifying
 * all elements that are within a maximum radius (a circle shape).
 */
MatrixXb circleMask(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid,
                    double origin[2], double radius){
    MatrixXd rgrid;  // double
    MatrixXb mask;   // boolean

    // squared euclidean distance (not taking sqrt, square the radius instead)
    rgrid = (xgrid.array() - origin[0]).pow(2) + (ygrid.array() - origin[1]).pow(2);
    // convert to boolean based on distance from origin vs radius of desired circle
    mask = rgrid.array() <= pow(radius, 2);
    return mask;
}

// Rotate X and Y grid matrices clockwise by given degrees, around an origin.
std::tuple<MatrixXd, MatrixXd> rotateGrids(double origin[2], MatrixXd xgrid, MatrixXd ygrid,
                                           double degrees){
    MatrixXd x_rot, y_rot;  // double

    // convert orientation angle to radians
    double theta = degrees * 3.14159265359/180;
    // rotate x and y grids around origin
    x_rot = origin[0] + cos(theta)*(xgrid.array()-origin[0]) - sin(theta)*(ygrid.array()-origin[1]);
    y_rot = origin[1] + sin(theta)*(xgrid.array()-origin[0]) - cos(theta)*(ygrid.array()-origin[1]);

    std::tuple<MatrixXd, MatrixXd> out = std::make_tuple(x_rot, y_rot);
    return out;
}

/* Take X and Y grid matrices and draw a rectangular boolean mask. Rectangle is defined by WxH dims
 * centred on an origin, and oriented in an arbitrary angle.
 */
MatrixXb rectMask(MatrixXd xgrid, MatrixXd ygrid, double origin[2], double orient,
                  double width, double height){
    MatrixXb mask; // boolean

    // rotate coordinates according to orientation of desired rectangle mask
    auto [xgrid_rot, ygrid_rot] = rotateGrids(origin, xgrid, ygrid, orient);
    // convert to boolean based on distances from origin on rotated x and y planes
    mask = ((xgrid_rot.array() - origin[0]).abs() <= width/2) * ((ygrid_rot.array() - origin[1]).abs() <= height/2);
    return mask;
}
