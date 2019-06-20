//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_OFFALPHA_H
#define RETINA_OFFALPHA_H

#include <random>

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"


class OffAlpha : public Cell {
protected:
    double tonic;

public:
    OffAlpha(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
            Eigen::VectorXd &yOnes, const double net_dt, const std::array<double, 2> cell_pos, std::mt19937 gen)
            :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "OffAlpha";
        // spatial properties
        diam = 8;  // 15
        rf_rad = 100;  // 200
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        tonic = 3.14159265359 * pow(rf_rad, 2);  // area of the receptive field (for sustained)
        onoff = false;
        sustained = !std::uniform_int_distribution<> (0, 1)(gen);  // randomly set cell to sustained or transient
        dtau = sustained ? 25 : 100;  // decay tau depends on sustained/transient status (true : false)
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
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>() * -1;
        return mask;
    }

    void stimulate(double strength, double angle) override {
        /* Since the stim representations are sparse, I don't want the default to be -ve,
         * so instead, a strength=0 passed from the Stim object will result in excitation.
         */
        // tonic stimulation only for sustained cells. Scale stim strength accordingly
        Vm += sustained*tonic*.02 + strength*(sustained ? .24 : 4);
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


#endif //RETINA_OFFALPHA_H
