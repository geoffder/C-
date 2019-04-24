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
              const double cell_pos[2], const double cell_diam, const double rf, const double cell_dtau)
              :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos, cell_diam, cell_dtau) {
        /* it has occurred to me that while the inputs to the Cell constructor made sense initially, I think
         * that I should change the way some of the parameters are set. For instance, each cell type will have a
         * given dtau and soma size, so those shouldn't have to be passed as a variable into the constructor when
         * net.populate() calls it. Instead, these derived cell classes will have defaults for those values that
         * they can set the inherited variables of Cell too.
         */
    }
};


#endif //RETINA_ONOFFDSGC_H
