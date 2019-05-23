//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_LOCALEDGEDETECTOR_H
#define RETINA_LOCALEDGEDETECTOR_H

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"


class LocalEdgeDetector : public Cell {
public:
    LocalEdgeDetector(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
                        Eigen::VectorXd &yOnes, const double net_dt, const double cell_pos[2])
                        :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "LocalEdgeDetector";
        // spatial properties
        diam = 15;
        rf_rad = 45;
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = false;
        onoff = true;
        dtau = 200;
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
        Vm += strength*4;
    }
};


#endif //RETINA_LOCALEDGEDETECTOR_H
