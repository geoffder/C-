//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_ONALPHA_H
#define RETINA_ONALPHA_H

#include <random>

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class OnAlpha : public Cell {
public:
    OnAlpha(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
            Eigen::VectorXd &yOnes, const double net_dt, const std::array<double, 2> cell_pos, std::mt19937 gen)
            :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "OnAlpha";
        // spatial properties
        diam = 8;  // 15
        rf_rad = 100;  // 200
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        onoff = false;
        sustained = !std::uniform_int_distribution<> (0, 1)(gen);  // randomly set cell to sustained or transient
        dtau = sustained ? 200 : 100;  // decay tau depends on sustained/transient status (true : false)
    }

    Eigen::MatrixXi buildRF(Eigen::VectorXd xgrid, Eigen::VectorXd ygrid, Eigen::VectorXd xOnes,
                            Eigen::VectorXd yOnes, std::array<double, 2> origin, double radius) {
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (
                    (xgrid.array() - origin[0]).square().matrix() * yOnes.transpose()
                     + xOnes * (ygrid.array() - origin[1]).square().matrix().transpose()
                 );
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>();
        return mask;
    }

    void stimulate(double strength, double angle) override {
        // modulate stimulus strength depending on whether cell is sustained or transient
        Vm += strength*(sustained ? .025 : .5);
    }

    // override base method to add in sustained/transient identifier
    std::string getParamStr() override {
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "sustained": )" << sustained << R"(, "diam": )" << diam;
        stream << R"(, "rf_rad": )" << rf_rad << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_ONALPHA_H
