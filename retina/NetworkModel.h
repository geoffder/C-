#ifndef RETINA_NETWORKMODEL_H
#define RETINA_NETWORKMODEL_H

#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>
#include <random>

#include "Cell.h"
#include "Stim.h"
#include "utils.h"
#include "eigen_types.h"


class NetworkModel {
private:
    int dims[2];
    Eigen::MatrixXd xgrid;
    Eigen::MatrixXd ygrid;
    int margin;
    int tstop;
    double dt;                            // timestep of network
    double t;
    int runs;

    std::vector<Cell> cells;
    std::vector<Stim> stims;

public:
    NetworkModel(const int net_dims[2], const int cell_margin, const int time_stop, const double delta){
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        auto [xgrid, ygrid] = gridMats(dims[0], dims[1]);
        margin = cell_margin;
        tstop = time_stop;
        dt = delta;                            // timestep of network
        t = 0;
        runs = 0;
    }

    void populate(const int spacing, double jitter){
        double theta, radius;
        double pos[2];
        Eigen::VectorXd xgridvec = VectorXd::LinSpaced((dims[1]-margin*2)/spacing, margin, dims[1]-margin);
        Eigen::VectorXd ygridvec = VectorXd::LinSpaced((dims[1]-margin*2)/spacing, margin, dims[1]-margin);

        // random number generation
        std::random_device rd;  // Obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        // uniform distribution for angle cell is offset from grid-point
        std::uniform_real_distribution<> Uniform(0, 1.0);
        // Normal distribution for distance cell is offset from grid-point
        std::normal_distribution<> Gaussian{0, 1};


        for (int i; i < xgridvec.size(); ++i){
            for (int j; j < ygridvec.size(); ++j){
                theta = 2 * 3.14159265359 * Uniform(gen);
                radius = Gaussian(gen) * jitter;
                pos[0] = xgridvec[i] + radius * cos(theta);
                pos[1] = ygridvec[j] + radius * sin(theta);
                // make a new cell and push it to the list
            }
        }
    }

    void clearCells(){
        cells.clear();
    }
    void addStim(const Stim &new_stim){
        stims.push_back(new_stim);
    }

    void clearStims(){
        stims.clear();
    }
};



#endif //RETINA_NETWORKMODEL_H
