#ifndef RETINA_NETWORKMODEL_H
#define RETINA_NETWORKMODEL_H
//
//#include <iostream>
//#include <Eigen/Dense>
//#include <tuple>
//#include <vector>

#include "Cell.h"
#include "Stim.h"
//#include "utils.h"
//#include "eigen_types.h"
#include "everything.h"
//using namespace utils;
//using namespace Eigen;
//typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

//std::tuple<MatrixXd, MatrixXd> gridMats(int nrows, int ncols);

class NetworkModel {
private:
    std::tuple<int, int> dims;
    Eigen::MatrixXd xgrid;
    Eigen::MatrixXd ygrid;
    int margin;
    int tstop;
    float dt;                            // timestep of network
    int t;
    int runs;

    std::vector<Cell> cells;
    std::vector<Stim> stims;
};


#endif //RETINA_NETWORKMODEL_H
