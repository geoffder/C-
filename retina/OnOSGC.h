//
// Created by geoff on 2019-05-10.
//

#ifndef RETINA_ONOSGC_H
#define RETINA_ONOSGC_H

#include <math.h>
#include <random>

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"


class OnOSGC : public Cell {
protected:
    double theta;      // preferred orientation
    double axis0;
    double axis1;
    std::array<double, 3> cardinals = {0, 90};  // constant

public:
    OnOSGC(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const double net_dt,
           const double cell_pos[2], std::mt19937 gen)
            :Cell(net_dims, xgrid, ygrid, net_dt, cell_pos) {
        type = "OnOSGC";
        // spatial properties
        diam = 15;  // of soma
        somaMask = circleMask(*net_xgrid, *net_ygrid, pos, diam/2);
        // Orientation-selective properties
        axis0 = 50;
        axis1 = 150;
        theta = rollPreferred(gen);  // choose a cardinal direction preference for this cell
        rfMask = buildRF(*net_xgrid, *net_ygrid, pos, axis0, axis1, theta);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        onoff = false;
        dtau = 250;
    }

    double rollPreferred(std::mt19937 gen) {
        // sample one cardinal direction from the array
        std::vector<double> choice;  // iterator to receive sample output
        std::sample(cardinals.begin(), cardinals.end(), std::back_inserter(choice), 1, gen);
        return choice[0];
    }

    // axis0 and axis1 are the full length of the minor and major axes of the ellipse (like diam, not rad)
    Eigen::MatrixXi buildRF(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double axis0,
                            double axis1, double theta) {
        Eigen::MatrixXd x, y;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        x = (xgrid.array() - origin[0])*cos(theta) + (ygrid.array() - origin[1])*sin(theta);
        y = (xgrid.array() - origin[0])*sin(theta) + (ygrid.array() - origin[1])*cos(theta);
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (((x.array()/axis0).square() + (y.array()/axis1).square()) <= 1).cast<int>();
        return mask;
    }

    // override base method to add in theta
    std::string getParamStr() override {
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "theta": )" << theta << R"(, "diam": )" << diam;
        stream << R"(, "axis0": )" << axis0 << R"(, "axis1": )" << axis1 << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_ONOSGC_H
