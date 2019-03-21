#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include "types.h"

using namespace eigen_types;
using namespace Eigen;
//typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

/* This project is for fiddling with Eigen and getting familiar
 * with its API. Broadcasting and creation of masks in specific.
 * For now I'll be accumulating a novice """cookbook""" in here
 * as I learn the ropes.
 */


/* This function takes the dimensions (y, x) and returns a tuple packed with
 * two matrices. The X and Y matrices serve a similar purpose to numpy's ogrid
 * function, which returns an arange column vector and an arange row vector.
 * These can be used for the creation of spatial masks.
 */
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
MatrixXb circleMask(const MatrixXd &xgrid, const MatrixXd &ygrid, const double origin[2], const double radius){
    MatrixXd rgrid;  // double
    MatrixXb mask;   // boolean

    // squared euclidean distance (not taking sqrt, square the radius instead)
    rgrid = (xgrid.array() - origin[0]).pow(2) + (ygrid.array() - origin[1]).pow(2);
    // convert to boolean based on distance from origin vs radius of desired circle
    mask = rgrid.array() <= pow(radius, 2);
    return mask;
}

// Rotate X and Y grid matrices clockwise by given degrees, around an origin.
std::tuple<MatrixXd, MatrixXd> rotateGrids(const double origin[2], const MatrixXd &xgrid, const MatrixXd &ygrid,
                                           const double degrees){
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
MatrixXb rectMask(const MatrixXd &xgrid, const MatrixXd &ygrid, const double origin[2], const double orient,
                  const double width, const double height){
    MatrixXb mask; // boolean

    // rotate coordinates according to orientation of desired rectangle mask
    auto [xgrid_rot, ygrid_rot] = rotateGrids(origin, xgrid, ygrid, orient);
    // convert to boolean based on distances from origin on rotated x and y planes
    mask = ((xgrid_rot.array() - origin[0]).abs() <= width/2) * ((ygrid_rot.array() - origin[1]).abs() <= height/2);
    return mask;
}

int test_ops() {
    MatrixXi m;  // dynamic sized matrix
    m = MatrixXi::Ones(10, 4)*5;  //  fill with (10 x 4) shape of ones
    VectorXi v(m.outerSize());  // vector as long as number of cols of m
    v << 0, 1, 2, 3;  // stream in values for vector
    std::cout << v << "\n";

    // broadcast vector over every row of the matrix
    m = m.array().rowwise() * v.transpose().array();
    m(1, 2) = 10;  // set value at index
    std::cout << m(1, 2) << "\n\n";
    std::cout << m << "\n\n";

    int nrows, ncols;
    nrows = 10, ncols = 10;
    // auto takes the returned tuple and automatically creates variables with
    // the given names of the proper types to unpack it.
    auto [xgrid, ygrid] = gridMats(nrows, ncols);

    double origin[2] = {5, 5};
    double radius = 2.2;
    MatrixXb circle = circleMask(xgrid, ygrid, origin, radius);
    std::cout << circle << "\n\n";

    double orient, width, height;
    orient = 45, width = 5, height = 2;
    MatrixXb rectangle = rectMask(xgrid, ygrid, origin, orient, width, height);
    std::cout << rectangle << "\n";
    return 0;
}


