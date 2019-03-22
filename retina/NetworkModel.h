#ifndef RETINA_NETWORKMODEL_H
#define RETINA_NETWORKMODEL_H

#include <iostream>
#include <tuple>
#include <vector>
#include <chrono>
#include <random>

#include <Eigen/Dense>

#include "Cell.h"
#include "Stim.h"
#include "utils.h"
#include "eigen_types.h"


class NetworkModel {
private:
    int dims[2];
    Eigen::MatrixXd xgrid;                  // coordinate range (along rows) grid used to calculate masks
    Eigen::MatrixXd ygrid;                  // coordinate range (down columns) grid used to calculate masks
    int margin;                             // cell free margin (allow stimulus to move in from outside)
    int tstop;                              // end time of a simulation run
    double dt;                              // timestep of network
    double t;                               // current time
    int runs;                               // number of runs completed in this Network simulation

    std::vector<Cell> cells;                // All Cell objects residing in this Network
    std::vector<double> cell_Xs;            // Centre X coordinates of all cells in the network
    std::vector<double> cell_Ys;            // Centre Y coordinates of all cells in the network
    std::vector<Stim> stims;                // All Stimuli objects added to the Network

public:
    NetworkModel(const int net_dims[2], const int cell_margin, const int time_stop, const double delta){
        // spatial
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        std::tie(xgrid, ygrid) = gridMats(dims[0], dims[1]);
        margin = cell_margin;
        // temporal
        tstop = time_stop;
        dt = delta;
        t = 0;
        runs = 0;
    }

    void populate(const int spacing, double jitter){
        double theta, radius;
        double pos[2];

        // Regularly spaced axes giving the default layout of cells in the network.
        Eigen::VectorXd xgridvec = Eigen::VectorXd::LinSpaced((dims[1]-margin*2)/spacing, margin, dims[1]-margin);
        Eigen::VectorXd ygridvec = Eigen::VectorXd::LinSpaced((dims[1]-margin*2)/spacing, margin, dims[1]-margin);

        // Random number generation
        std::random_device rd;  // Obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        // uniform distribution for angle cell is offset from grid-point
        std::uniform_real_distribution<> Uniform(0, 1.0);
        // Normal distribution for distance cell is offset from grid-point
        std::normal_distribution<> Gaussian{0, 1};

        for (int i = 0; i < xgridvec.size(); ++i){
            for (int j = 0; j < ygridvec.size(); ++j){

                theta = 2 * 3.14159265359 * Uniform(gen);
                radius = Gaussian(gen) * jitter;

                pos[0] = xgridvec[i] + radius * cos(theta);
                pos[1] = ygridvec[j] + radius * sin(theta);

                // make a new cell and push it to the list (default cells for now, need to create
                // different cell types to pick from. Sub-classes of cell with their own default parameters.
                // emplace_back is push_back, but it constructs a new object of the vector's type
                cells.emplace_back(dims, xgrid, ygrid, dt, pos, double(15), double(15), double(.5));
                cell_Xs.push_back(pos[0]), cell_Ys.push_back(pos[1]);
            }
        }
    }

    std::vector<Cell> getCells(){
        return cells;
    }

    // Cell Coordinates organized into a Nx2 matrix (for output to file)
    Eigen::Matrix<double, Eigen::Dynamic, 2> getCellXYs(){
        Eigen::MatrixXd cell_XYs = Eigen::MatrixXd::Ones(cell_Xs.size(), 2);
        for(std::size_t i = 0; i < cell_Xs.size(); ++i){
            cell_XYs(i, 0) = cell_Xs[i];
            cell_XYs(i, 1) = cell_Ys[i];
        }
        return cell_XYs;
    }

    // Construct table with columns of all cell recordings (for output to file)
    Eigen::MatrixXd getRecTable(){
        // work on a better way using matrix concatenation in Eigen! (Like below example)
        // MatrixXd C(A.rows(), A.cols()+B.cols());
        // C << A, B;
        // need to get the std::vectors in to Eigen Matrices/Vectors first of course
        int recLen = cells[0].getRecs()[0].size() * cells[0].getRecs().size();
        Eigen::MatrixXd all_recs = Eigen::MatrixXd::Ones(recLen, cells.size());
        for(std::size_t c = 0; c < cells.size(); ++c){
            auto recs = cells[c].getRecs();
            for(std::size_t i = 0; i < recs.size(); ++i){
                for(std::size_t j = 0; j < recs[i].size(); ++j){
                    all_recs(i*recs[i].size()+j, c) = recs[i][j];
                }
            }
        }
        return all_recs;
    }

    void clearCells(){
        cells.clear();
    }

    // Generating new Stimuli within NetworkModel for ease of access to it's variables.
    void newStim(const double start_pos[2], const int time_on, const int time_off, const double velocity,
                 const double direction, const double orientation, const double amplitude, const double change,
                 const std::string &type, const double radius=50, const double width=50, const double length=100){

        stims.emplace_back(dims, xgrid, ygrid, dt, start_pos, time_on, time_off, velocity, direction, orientation,
                amplitude, change);

        if (type == "bar"){
            stims[stims.size()-1].setBar(width, length);
        } else if (type == "circle"){
            stims[stims.size()-1].setCircle(radius);
        }
    }

    std::vector<Stim> getStims(){
        return stims;
    }

    void clearStims(){
        stims.clear();
    }

    void step(){
        double strength;
        Eigen::SparseMatrix<int> * sparseRF_ref;

        for(auto& stim : stims){
            stim.move();
            for(auto& cell : cells){
                sparseRF_ref = cell.getSparseRFref();
                strength = stim.check(sparseRF_ref);
                cell.excite(strength);
            }
        }

        for(auto& cell : cells){
            cell.decay();
        }
        t += dt;
    }

    void run(){
        auto run_start = Clock::now();
        t = 0;
        for(int i = 0; i*dt < tstop; ++i){
            step();
        }
        for(auto& cell : cells){
            cell.storeRec();
            cell.setVm(0);
        }
        ++runs;
        auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - run_start).count();
        std::cout << run_time << " milliseconds run time\n";
    }

    // Add up spatial masks of all cells and their receptive fields for display (for output to file)
    Eigen::MatrixXd cellMatrix(){
        Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(dims[0], dims[1]);
        for(auto& cell : cells){
            mat += cell.getSoma().cast<double>();
            mat += cell.getRF().cast<double>()*.2;
        }
        return mat;
    }
};



#endif //RETINA_NETWORKMODEL_H