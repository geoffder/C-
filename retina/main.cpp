//#include <iostream>
//#include <Eigen/Dense>
//#include <tuple>
//#include <vector>
//
//#include "NetworkModel.h"
//#include "Cell.h"
//#include "Stim.h"
//#include "eigen_types.h"
//#include "NetworkModel.h"

//#include "utils.h"
#include "everything.h"
using namespace utils;
using namespace Eigen;
//typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

//std::tuple<MatrixXd, MatrixXd> gridMats(int nrows, int ncols);
//MatrixXb circleMask(MatrixXd &xgrid, MatrixXd &ygrid, double origin[2], double radius);

//std::tuple<MatrixXd, MatrixXd> rotateGrids(double origin[2], MatrixXd &xgrid, MatrixXd &ygrid, double degrees);
//MatrixXb rectMask(MatrixXd &xgrid, MatrixXd &ygrid, double origin[2], double orient, double width, double height);

//int test_ops();

int main() {
    std::cout << "Hello, World!" << std::endl;

    int dims[2] = {10, 10};
    int dt = 1;
    double pos[2] = {5, 5};
    double diam = 5;
    double rf = 7;
    float dtau = 1;

//    MatrixXi test = MatrixXi::Ones(10, 10);
//    std::cout << test;
//    auto [xgrid, ygrid] = gridMats(dims[0], dims[1]);
//    std::cout << xgrid << "\n";

//    MatrixXb circ = xgrid.array() > 4;
//    std::cout << circ;
//    MatrixXb circ = circleMask(xgrid, ygrid, pos, diam/2);
//    Cell a(dims, xgrid, ygrid, dt, pos, diam, rf, dtau);

//    a.setVm(10);
//    std::cout << "voltage: " << a.getVm() << "\n";
//
//    for(int i = 0; i < 10; ++i){
//      a.decay();
//    }
//
//    std::vector<float> recording = a.getRec();
//
//    for(int i = 0; i < 10; ++i){
//        std::cout << recording[i] << " ";
//    }
//    std::cout << "\n";

//    test_ops();

    return 0;
}