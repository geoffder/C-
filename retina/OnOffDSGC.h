//
// Created by geoff on 2019-04-23.
//

#ifndef RETINA_ONOFFDSGC_H
#define RETINA_ONOFFDSGC_H

#include <math.h>
#include <random>

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class OnOffDSGC : protected Cell {
protected:
    double prefInhib;  // minimum inhibitory input mod (when stim angle matches preferred angle)
    double nullInhib;  // maximum inhibitory input mod (when stim angle matches null angle)
    double theta;      // preferred angle
    std::array<double, 4> cardinals= {0, 90, 180, 270};  // constant

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
        // active / synaptic properties
        sustained = false;
        dtau = .5;
        // Direction-selective properties
        type = "OnOffDSGC";
        prefInhib = 0;  // DSGC specific constant (inhibition in preferred direction)
        nullInhib = 1.5;  // DSGC specific constant (inhibition in null direction)
        theta = rollPreferred();  // choose a cardinal direction preference for this cell
    }

    double rollPreferred() {
        // Random number generation (consider passing the cell constructor the generator from populate)
        std::random_device rd;  // Obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        // sample one cardinal direction from the array
        std::vector<double> choice;  // iterator to receive sample output
        std::sample(cardinals.begin(), cardinals.end(), std::back_inserter(choice), 1, gen);
        return choice[0];
    }

    Eigen::MatrixXi buildRF(Eigen::MatrixXd xgrid, Eigen::MatrixXd ygrid, double origin[2], double radius){
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi mask;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (xgrid.array() - origin[0]).square() + (ygrid.array() - origin[1]).square();
        // convert to boolean based on distance from origin vs radius of desired circle
        mask = (rgrid.array() <= pow(radius, 2)).cast<int>();
        return mask;
    }

    void stimulate(double strength, double angle) override {
        excite(strength);
        inhibit(strength, angle);
    }

    void excite(double strength) {
        Vm += strength;
    }

    void inhibit(double strength, double angle) {
        // can take Vm below zero, but decay step will clip to zero
        Vm -= strength * (prefInhib + (nullInhib - prefInhib) * tuning(angle));
    }

    // sigmoid scaling from ~ 0 -> 1 (strength mod), over 0 -> 180 in degrees off preferred
    double tuning(double angle) {
        return 1.0 - .98/(1.0 + std::exp((angle - 91.0)/25.0));
    }

    // override base method to add in theta
    std::string getParamStr() override {
        std::stringstream stream;
        // JSON formatting using raw string literals
        stream << R"({"type": ")" << type << R"(", "theta": )" << theta << type << R"(", "diam": )" << diam;
        stream << R"(, "rf_rad": )" << rf_rad << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_ONOFFDSGC_H
