//
// Created by geoff on 2019-04-23.
//

#ifndef RETINA_ONOFFDSGC_H
#define RETINA_ONOFFDSGC_H

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class OnOffDSGC : protected Cell {

public:
    OnOffDSGC(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
              const double cell_pos[2])
              :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos) {
        /* it has occurred to me that while the inputs to the Cell constructor made sense initially, I think
         * that I should change the way some of the parameters are set. For instance, each cell type will have a
         * given dtau and soma size, so those shouldn't have to be passed as a variable into the constructor when
         * net.populate() calls it. Instead, these derived cell classes will have defaults for those values that
         * they can set the inherited variables of Cell too.
         */
        // spatial properties
        diam = 15;
        rf_rad = 60;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = buildRF(*net_xgrid, *net_ygrid, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert to sparse, everything below 1 is zeroed (not anymore?)
        // active properties
        dtau = .5;
    }

    Eigen::MatrixXi buildRF(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius){
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi mask;   // boolean

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (xgrid.array() - origin[0]).square() + (ygrid.array() - origin[1]).square();
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>();
        return mask;
    }

    // override base class method
    // ON-OFF, so needs to respond only to onsets. DS so needs to take in angle.
    // thus, need to add theta as parameter to ALL excite methods (called from Stim) even if not used by the cell.
    /* New idea from talking to Ben:
     * When I draw the mask for each stimuli on each timestep, also create a DIFFERENCE mask (with sparse version)
     * that describes the "pixel" values that have changed (negative or positive) while all non-changed pixels are
     * zero.
     */
    void excite(double strength, double theta){
        Vm += strength;
    }
};


#endif //RETINA_ONOFFDSGC_H
