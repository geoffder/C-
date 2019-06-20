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
    //Eigen::MatrixXd * net_xgrid;                 // pointer to network X range grid used for generation of masks
    //Eigen:: MatrixXd * net_ygrid;                // pointer to network Y range grid used for generation of masks
    Eigen::VectorXd * net_xvec;                 // pointer to network X range grid used for generation of masks
    Eigen::VectorXd * net_yvec;                // pointer to network Y range grid used for generation of masks
    Eigen::VectorXd * net_xOnes;                 // pointer to network X range grid used for generation of masks
    Eigen::VectorXd * net_yOnes;
    double dt;                                   // timestep of network model
    // cell and spatial properties
    std::string type = "base";
    std::array<double, 2> pos;                               // centre coordinates (constant)
    double diam;                                 // soma diameter
    double rf_rad;                               // receptive field radius
    Eigen::MatrixXi somaMask;                    // mask defining cell body
    Eigen::MatrixXi rfMask;                      // mask defining receptive field
    Eigen::SparseMatrix<int> rfMask_sparse;      // sparse representation of the receptive field (fast computation)
    // active properties
    bool sustained;                              // whether cell is sustained (otherwise transient)
    bool onoff;                                  // whether cell is an OnOff cell
    double Vm;                                   // "membrane" state
    double dtau;                                 // decay tau
    std::vector<double> rec;                     // activity recording

public:
    // default constructor (should never be used, just delete I think)
    Cell() {
        Vm = 0;
    }

    // network and generic cell properties only, used by derived classes
    Cell(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
            Eigen::VectorXd &yOnes,const double net_dt, const std::array<double, 2> cell_pos) {
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xvec = &xgrid;  // point to the xgrid of the Network this cell belongs to (memory efficiency)
        net_yvec = &ygrid;  // point to the ygrid of the Network this cell belongs to
        net_xOnes = &xOnes;
        net_yOnes = &yOnes;
        dt = net_dt;
        // spatial properties (
        pos[0] = cell_pos[0], pos[1] = cell_pos[1];
        // active properties
        Vm = 0;
    }

    Eigen::MatrixXi getSoma() {
        return somaMask;
    }

    Eigen::MatrixXi getRF() {
        return rfMask;
    }

    Eigen::SparseMatrix<int>* getSparseRFref() {
        return &rfMask_sparse;
    }

    bool isSustained() {
        return sustained;
    }

    bool isOnOff() {
        return onoff;
    }

    double getVm() {
        return Vm;
    }

    void setVm(double newVm) {
        Vm = newVm;
    }

    std::vector<double> getRec() {
        return rec;
    }

    void clearRec() {
        rec.clear();
    }

    virtual void stimulate(double strength, double angle) {
        // angle used for only some cell types
        Vm += strength;
    }

    virtual void decay() {
        double delta = Vm * (1 - exp(-dt/dtau));
        Vm = std::fmax(double (0), Vm - delta);
        rec.push_back(Vm);
    }

    virtual std::string getParamStr() {
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "diam": )" << diam << R"(, "rf_rad": )" << rf_rad;
        stream << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_CELL_H
