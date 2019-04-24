#ifndef RETINA_CELL_H
#define RETINA_CELL_H

#include <iostream>
#include <tuple>
#include <vector>
#include <chrono>
#include <sstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "utils.h"
#include "type_defs.h"


class Cell {
protected:
    // network properties
    int dims[2];                                 // dimensions of network model this cell belongs to
    Eigen::MatrixXd * net_xgrid;                 // pointer to network X range grid used for generation of masks
    Eigen:: MatrixXd * net_ygrid;                // pointer to network Y range grid used for generation of masks
    double dt;                                   // timestep of network model
    // cell and spatial properties
    std::string type = "base";
    double pos[2];                               // centre coordinates (constant)
    double diam;                                 // soma diameter
    double rf_rad;                               // receptive field radius
    Eigen::MatrixXi somaMask;                    // mask defining cell body
    Eigen::MatrixXi rfMask;                      // mask defining receptive field
    Eigen::SparseMatrix<int> rfMask_sparse;      // sparse representation of the receptive field (fast computation)
    // active properties
    double Vm;                                   // "membrane" state
    double dtau;                                 // decay tau
    std::vector<double> rec;                     // activity recording

public:
    // default constructor
    Cell() {
        Vm = 0;
    }
    // network and generic cell properties only, used by derived classes
    Cell(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
         const double cell_pos[2], const double cell_diam) {
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xgrid = &xgrid;  // point to the xgrid of the Network this cell belongs to (memory efficiency)
        net_ygrid = &ygrid;  // point to the ygrid of the Network this cell belongs to
        dt = net_dt;
        // spatial properties
        pos[0] = cell_pos[0], pos[1] = cell_pos[1];
        diam = cell_diam;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        // active properties
        Vm = 0;
    }
    Cell(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
            const double cell_pos[2], const double cell_diam, const double rf, const double cell_dtau){
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xgrid = &xgrid;  // point to the xgrid of the Network this cell belongs to (memory efficiency)
        net_ygrid = &ygrid;  // point to the ygrid of the Network this cell belongs to
        dt = net_dt;
        // spatial properties
        pos[0] = cell_pos[0], pos[1] = cell_pos[1];
        diam = cell_diam;
        rf_rad = rf;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = circleMask(*net_xgrid, *net_ygrid, pos, rf);
        rfMask_sparse = rfMask.sparseView();  // convert to sparse, everything below 1 is zeroed (not anymore?)
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

    void clearRec(){
        rec.clear();
    }

    void excite(double strength){
        Vm += strength;
    }

    void decay(){
        double delta = Vm * (1 - exp(-dt/dtau));
        Vm = std::fmax(double (0), Vm - delta);
        rec.push_back(Vm);
    }

    std::string getParamStr(){
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "diam": )" << diam << R"(, "rf_rad": )" << rf_rad;
        stream << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_CELL_H
