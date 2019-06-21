//
// Created by geoff on 2019-04-23.
//

#ifndef RETINA_ONOFFDSGC_H
#define RETINA_ONOFFDSGC_H

#include <array>
#include <math.h>
#include <random>

#include "Cell.h"
#include "utils.h"
#include "type_defs.h"

class OnOffDSGC : public Cell {
protected:
    double prefInhib;  // minimum inhibitory input mod (when stim angle matches preferred angle)
    double nullInhib;  // maximum inhibitory input mod (when stim angle matches null angle)
    double theta;      // preferred angle
    std::array<double, 4> cardinals = {0, 90, 180, 270};  // constant

public:
    OnOffDSGC(const std::array<int, 2> net_dims, Eigen::VectorXd &xgrid, Eigen::VectorXd &ygrid, Eigen::VectorXd &xOnes,
                Eigen::VectorXd &yOnes, const double net_dt, const std::array<double, 2> cell_pos, std::mt19937 gen)
                :Cell(net_dims, xgrid, ygrid, xOnes, yOnes, net_dt, cell_pos) {
        type = "OnOffDSGC";
        // spatial properties
        diam = 8; // 15
        centre_rad = 50; // 100
        surround_rad = 100;
        somaMask = circleMask(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, diam/2);
        std::tie(rfCentre, rfSurround) = buildRF(*net_xvec, *net_yvec, *net_xOnes, *net_yOnes, pos, centre_rad, surround_rad);
        rfCentre_sparse = rfCentre.sparseView();  // convert from dense matrix to sparse
        rfSurround_sparse = rfSurround.sparseView();  // convert from dense matrix to sparse
        // active / synaptic properties
        sustained = false;
        onoff = true;
        dtau = 100;
        // Direction-selective properties
        prefInhib = 0;  // DSGC specific constant (inhibition in preferred direction)
        nullInhib = 1.5;  // DSGC specific constant (inhibition in null direction)
        theta = rollPreferred(gen);  // choose a cardinal direction preference for this cell
    }

    double rollPreferred(std::mt19937 gen) {
        // sample one cardinal direction from the array
        std::vector<double> choice;  // iterator to receive sample output
        std::sample(cardinals.begin(), cardinals.end(), std::back_inserter(choice), 1, gen);
        return choice[0];
    }

    rfPair buildRF(Eigen::VectorXd xgrid, Eigen::VectorXd ygrid, Eigen::VectorXd xOnes,
                            Eigen::VectorXd yOnes, std::array<double, 2> origin, double radius, double surradius) {
        Eigen::MatrixXd rgrid;  // double
        Eigen::MatrixXi centre, surround;   // integer

        // squared euclidean distance (not taking sqrt, square the radius instead)
        rgrid = (
                (xgrid.array() - origin[0]).square().matrix() * yOnes.transpose()
                + xOnes * (ygrid.array() - origin[1]).square().matrix().transpose()
        );
        // convert to boolean based on distance from origin vs radius of desired circle
        centre = (rgrid.array() <= pow(radius, 2)).cast<int>();
        surround = (pow(radius, 2) <= rgrid.array()).cast<int>()
                   * (rgrid.array() <= pow(radius+surradius, 2)).cast<int>();

        return std::make_tuple(centre, surround);
    }

    void check(Stim &stim) override {
        // use stimulus itself if sustained, and delta of stimulus if transient
        Eigen::SparseMatrix<int> stim_mask = sustained ? stim.getSparseMask() : stim.getSparseDelta();

        Eigen::SparseMatrix<int> centre_overlap, surround_overlap;
        centre_overlap = stim_mask.cwiseProduct(rfCentre_sparse);
        surround_overlap = stim_mask.cwiseProduct(rfSurround_sparse);

        double centre_strength = centre_overlap.sum() * stim.getAmp();
        double surround_strength =  surround_overlap.sum() * stim.getAmp();

        double difference =  std::abs(theta-stim.getTheta());
        difference = difference < 180 ? difference : std::abs(difference - 360);

        Vm += centre_strength - surround_strength;
        inhibit(centre_strength, difference);
    }

    void stimulate(double strength, double angle) override {
        strength *= 1;
        double difference =  std::abs(theta-angle);
        if (difference > 180) {
            difference = std::abs(difference - 360);
        }
        excite(strength);
        inhibit(strength, difference);
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
        stream << R"({"type": ")" << type << R"(", "theta": )" << theta << R"(, "diam": )" << diam;
        stream << R"(, "centre_rad": )" << centre_rad << R"(, "surround_rad": )" << surround_rad << R"(, "dtau": )" << dtau << "}";
        std::string params = stream.str();
        return params;
    }
};


#endif //RETINA_ONOFFDSGC_H
