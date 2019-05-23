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
    OnOSGC(const int net_dims[2], Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
            Eigen::VectorXd &yOnes, const double net_dt, const double cell_pos[2], std::mt19937 gen)
            :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "OnOSGC";
        // spatial properties
        diam = 15;  // of soma
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        // Orientation-selective properties
        axis0 = 50;
        axis1 = 150;
        theta = rollPreferred(gen);  // choose a cardinal direction preference for this cell
        rfMask = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, axis0, axis1, theta);
        rfMask_sparse = rfMask.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = true;
        onoff = false;
        dtau = 200;
    }

    double rollPreferred(std::mt19937 gen) {
        // sample one cardinal direction from the array
        std::vector<double> choice;  // iterator to receive sample output
        std::sample(cardinals.begin(), cardinals.end(), std::back_inserter(choice), 1, gen);
        return choice[0];
    }

    // axis0 and axis1 are the full length of the minor and major axes of the ellipse (like diam, not rad)
    Eigen::MatrixXi buildRF(Eigen::VectorXd xgrid, Eigen::VectorXd ygrid, Eigen::VectorXd xOnes,
                            Eigen::VectorXd yOnes, double origin[2], double axis0, double axis1, double theta) {
        Eigen::MatrixXd x, y;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        xgrid = xgrid.array() - origin[0];
        ygrid = ygrid.array() - origin[1];
        x = xgrid*cos(theta) + ygrid*sin(theta);
        y = xgrid*sin(theta) + ygrid*cos(theta);

        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (
                    (
                            (x.array()/axis0).square().matrix()* yOnes.transpose()
                             + xOnes * (y.array()/axis1).square().matrix().transpose()
                    ).array() <= 1
                ).cast<int>();
        return mask;
    }

    // override base method to add in theta
    std::string getParamStr() override {
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "theta": )" << theta << R"(, "diam": )" << diam;
        stream << R"(, "rf_ax0": )" << axis0 << R"(, "rf_ax1": )" << axis1 << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }

    void stimulate(double strength, double angle) override {
        Vm += strength*.1;
    }
};


#endif //RETINA_ONOSGC_H
