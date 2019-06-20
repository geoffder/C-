//
// Created by geoff on 2019-05-03.
//

#ifndef RETINA_BASICCELL_H
#define RETINA_BASICCELL_H

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class BasicCell : public Cell {
public:
    BasicCell(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
                      Eigen::VectorXd &yOnes, const double net_dt, const std::array<double, 2> cell_pos)
            :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "Basic";
        // spatial properties
        diam = 8;  // 15
        rf_rad = 100;
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        onoff = false;
        dtau = 10;
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
};


#endif //RETINA_BASICCELL_H
