//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_OFFALPHA_H
#define RETINA_OFFALPHA_H

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"


class OffAlpha : public Cell {
protected:
    double tonic;

public:
    OffAlpha(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
            const double cell_pos[2])
            :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos) {
        type = "OffAlpha";
        // spatial properties
        diam = 15;
        rf_rad = 200;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = buildRF(*net_xgrid, *net_ygrid, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        tonic = 3.14159265359 * pow(rf_rad, 2);  // area of the receptive field
        sustained = true;
        onoff = false;
        dtau = 10;
    }

    Eigen::MatrixXi buildRF(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius) {
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (xgrid.array() - origin[0]).square() + (ygrid.array() - origin[1]).square();
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>() * -1;  // inhibited by positive stimuli
        return mask;
    }

    void stimulate(double strength, double angle) override {
        /* Since the stim representations are sparse, I don't want the default to be -ve,
         * so instead, a strength=0 passed from the Stim object will result in excitation.
         */
        Vm += tonic*.2 + strength;
    }

};


#endif //RETINA_OFFALPHA_H
