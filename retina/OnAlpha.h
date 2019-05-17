//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_ONALPHA_H
#define RETINA_ONALPHA_H


#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class OnAlpha : public Cell {
public:
    OnAlpha(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
            Eigen::VectorXd &yOnes, const double net_dt, const double cell_pos[2])
            :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "OnAlpha";
        // spatial properties
        diam = 15;
        rf_rad = 200;
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        onoff = false;
        dtau = 250;
    }

    Eigen::MatrixXi buildRF(Eigen::VectorXd xgrid, Eigen::VectorXd ygrid, Eigen::VectorXd xOnes,
                            Eigen::VectorXd yOnes, double origin[2], double radius) {
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
        Vm += strength*.025;
    }
};


#endif //RETINA_ONALPHA_H
