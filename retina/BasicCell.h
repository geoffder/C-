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
    BasicCell(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
              const double cell_pos[2])
            :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos) {
        // spatial properties
        diam = 15;
        rf_rad = 60;
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        rfMask = circleMask(*net_xgrid, *net_ygrid, pos, rf_rad);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        dtau = 10;
    }
};


#endif //RETINA_BASICCELL_H
