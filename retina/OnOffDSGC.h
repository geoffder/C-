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
              :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos, cell_diam) {


    }
};


#endif //RETINA_ONOFFDSGC_H
