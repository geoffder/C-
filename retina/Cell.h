#ifndef RETINA_CELL_H
#define RETINA_CELL_H

//#include <iostream>
//#include <Eigen/Dense>
//#include <tuple>
//#include <vector>
//
//#include "utils.h"
//#include "eigen_types.h"
#include "everything.h"
//using namespace utils;
using namespace Eigen;
//typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

//std::tuple<MatrixXd, MatrixXd> gridMats(int nrows, int ncols);
//MatrixXb circleMask(MatrixXd &xgrid, MatrixXd &ygrid, double origin[2], double radius);

class Cell {
private:
    int dims[2];                          // dimensions of network model this cell belongs to
    MatrixXd net_xgrid;                   // X range grid used for generation of masks
    MatrixXd net_ygrid;                   // Y range grid used for generation of masks
    float dt;                             // timestep of network model
    double pos[2];                        // centre coordinates (constant)
    double diam;                          // soma diameter
    double rf_rad;                        // receptive field radius
    MatrixXb somaMask;                    // mask defining cell body
    MatrixXb rfMask;                      // mask defining receptive field
    float Vm;                             // "membrane" state
    float dtau;                           // decay tau
    std::vector<float> rec;               // activity recording
    std::vector<std::vector<float>> recs; // collection of recordings (each trial)

public:
    Cell(const int net_dims[2], const MatrixXd &xgrid, const MatrixXd &ygrid, const int net_dt,
            const double cell_pos[2], const double cell_diam, const double rf, const float cell_dtau){
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xgrid = xgrid;
        net_ygrid = ygrid;
        dt = net_dt;
        // spatial properties
        pos[0] = cell_pos[0], pos[1] = cell_pos[1];
        diam = cell_diam;
        rf_rad = rf;
//        somaMask = utils::circleMask(net_xgrid, net_ygrid, pos, diam/2);
//        rfMask = utils::circleMask(net_xgrid, net_ygrid, pos, rf);
        // active properties
        Vm = 0;
        dtau = cell_dtau;
    }

    MatrixXb getSoma(){
        return somaMask;
    }

    MatrixXb getRF(){
        return rfMask;
    }

    float getVm(){
        return Vm;
    }

    void setVm(float newVm){
        Vm = newVm;
    }

    std::vector<float> getRec(){
        return rec;
    }

    void storeRec(){
        recs.push_back(rec);
        rec.clear();
    }

    std::vector<std::vector<float>> getRecs(){
        return recs;
    }

    void excite(float strength){
        Vm += strength;
    }

    void decay(){
        float delta = Vm * (1 - exp(-dt/dtau));
        Vm = std::fmax(float (0), Vm - delta);
        rec.push_back(Vm);
    }
};


#endif //RETINA_CELL_H
