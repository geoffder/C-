#ifndef RETINA_CELL_H
#define RETINA_CELL_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tuple>
#include <vector>

#include "utils.h"
#include "eigen_types.h"


class Cell {
private:
    int dims[2];                                 // dimensions of network model this cell belongs to
    Eigen::MatrixXd * net_xgrid;                 // pointer to network X range grid used for generation of masks
    Eigen:: MatrixXd * net_ygrid;                // pointer to network Y range grid used for generation of masks
    double dt;                                   // timestep of network model
    double pos[2];                               // centre coordinates (constant)
    double diam;                                 // soma diameter
    double rf_rad;                               // receptive field radius
    Eigen::MatrixXi somaMask;                    // mask defining cell body
    Eigen::MatrixXi rfMask;                      // mask defining receptive field
    Eigen::SparseMatrix<int> rfMask_sparse;
    double Vm;                                   // "membrane" state
    double dtau;                                 // decay tau
    std::vector<double> rec;                     // activity recording
    std::vector<std::vector<double>> recs;       // collection of recordings (each trial)

public:
    Cell(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
            const double cell_pos[2], const double cell_diam, const double rf, const double cell_dtau){
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xgrid = &xgrid;
        net_ygrid = &ygrid;
        dt = net_dt;
        // spatial properties
        pos[0] = cell_pos[0], pos[1] = cell_pos[1];
        diam = cell_diam;
        rf_rad = rf;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = circleMask(*net_xgrid, *net_ygrid, pos, rf);
        rfMask_sparse = rfMask.sparseView(1);
        // active properties
        Vm = 0;
        dtau = cell_dtau;
    }

    Eigen::MatrixXi getSoma(){
        return somaMask;
    }

    Eigen::MatrixXi getRF(){
        return rfMask;
    }

    Eigen::MatrixXi* getRFref(){
        return &rfMask;
    }

    Eigen::SparseMatrix<int>* getSparseRFref(){
        return &rfMask_sparse;
    }

    double getVm(){
        return Vm;
    }

    void setVm(double newVm){
        Vm = newVm;
    }

    std::vector<double> getRec(){
        return rec;
    }

    void storeRec(){
        recs.push_back(rec);
        rec.clear();
    }

    std::vector<std::vector<double>> getRecs(){
        return recs;
    }

    void excite(double strength){
        Vm += strength;
    }

    void decay(){
        double delta = Vm * (1 - exp(-dt/dtau));
        Vm = std::fmax(double (0), Vm - delta);
        rec.push_back(Vm);
    }
};


#endif //RETINA_CELL_H
