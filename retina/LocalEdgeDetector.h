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
    LocalEdgeDetector(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
                      const double cell_pos[2])
                      :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos) {
        type = "LocalEdgeDetector";
        // spatial properties
        diam = 15;
        rf_rad = 45;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = buildRF(*net_xgrid, *net_ygrid, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        onoff = true;
        dtau = 250;
    }

    Eigen::MatrixXi buildRF(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius) {
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (xgrid.array() - origin[0]).square() + (ygrid.array() - origin[1]).square();
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>();
        return mask;
    }
};


#endif //RETINA_LOCALEDGEDETECTOR_H
