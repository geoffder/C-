#ifndef RETINA_STIM_H
#define RETINA_STIM_H

#include <iostream>
#include <tuple>
#include <vector>
#include <chrono>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "utils.h"
#include "type_defs.h"


class Stim {
private:
    int dims[2];                             // dimensions of network model this cell belongs to
    Eigen::MatrixXd * net_xgrid;             // pointer to network X range grid used for generation of masks
    Eigen::MatrixXd * net_ygrid;             // pointer to network Y range grid used for generation of masks
    double dt;                               // timestep of network model
    double pos[2];                           // centre coordinates
    int tOn;                                 // time stimulus appears
    int tOff;                                // time stimulus turns off
    double vel;                              // velocity
    double theta, theta_rad;                 // angle of movement
    double orient, orient_rad;               // angle of orientation
    double amp;                              // intensity of stimulus
    double dAmp;                             // rate and direction of change in intensity of stimulus
    std::string type;                        // type of stimulus tag/label (e.g. circle, bar)
    double radius;                           // stimulus radius (circle type parameter)
    double length;                           // length (bar type parameter)
    double width;                            // width (bar type parameter)
    Eigen::MatrixXi mask;                    // mask defining this stimulus
    Eigen::SparseMatrix<int> mask_sparse;    // sparse representation of the stimulus mask (fast computation)
    std::vector<double> xPosRec;             // stored position from each timestep
    std::vector<double> yPosRec;             // stored position from each timestep
    std::vector<double> ampRec;              // stored amplitude from each timestep
    std::vector<double> orientRec;           // stored angle of orientation from each timestep

public:
    Stim(const int net_dims[2], Eigen::MatrixXd &xgrid, Eigen::MatrixXd &ygrid, const int net_dt,
            const double start_pos[2], const int time_on, const int time_off, const double velocity,
            const double direction, const double orientation, const double amplitude, const double change){
        // network properties
        dims[0] = net_dims[0], dims[1] = net_dims[1];
        net_xgrid = &xgrid;
        net_ygrid = &ygrid;
        dt = net_dt;
        // general stim properties
        pos[0] = start_pos[0], pos[1] = start_pos[1];
        tOn = time_on;
        tOff = time_off;
        // movement
        vel = velocity;
        theta = direction;
        theta_rad = theta*3.14159265359/180;
        // spatial orientation (only matters for radially asymmetric objects)
        orient = orientation;
        orient_rad = orient*3.14159265359/180;
        // "contrast"
        amp = amplitude;
        dAmp = change;
    }

    void setCircle(const double rad){
        type = "circle";
        radius = rad;
    }

    void setBar(const double wid, const double len){
        type = "bar";
        width = wid;
        length = len;
    }

    void setOrientation(const double angle){
        orient = angle;
    }

    void setVelocity(const double velocity){
       vel = velocity;
    }

    Eigen::MatrixXi getMask(){
        return mask;
    }

    void drawMask(){
        if (type == "bar") {
            mask = rectMask(*net_xgrid, *net_ygrid, pos, orient, width, length);
        } else if (type == "circle") {
            mask = circleMask(*net_xgrid, *net_ygrid, pos, radius);
        }
    }

    // update centre coordinate of stimulus (if moving), then redraw mask
    void move(){
        if (vel != 0) {
            pos[0] += vel / dt * cos(theta_rad);
            pos[1] += vel / dt * sin(theta_rad);
            drawMask();
            mask_sparse = mask.sparseView();  // mask.sparseView(1);
        }
        amp += dAmp;
        // record stimulus characteristics that are subject to change for movie reconstruction
        xPosRec.push_back(pos[0]);
        yPosRec.push_back(pos[1]);
        ampRec.push_back(amp);
        orientRec.push_back(orient);
    }

    double checkOld(Eigen::MatrixXi *rfMask){
        // time comparison of mask multiplication using sparse and dense matrices
        // sparse matrices (convert first from dense)
        Eigen::SparseMatrix<int> rfMaskSparse = (*rfMask).sparseView();
        Eigen::SparseMatrix<int> maskSparse = mask.sparseView();
        auto t1 = Clock::now();
        Eigen::SparseMatrix<int> sparse_overlap = maskSparse.cwiseProduct(rfMaskSparse);
        double sparse_sum = sparse_overlap.nonZeros();
        double sparse_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - t1).count();
        // old way with dense matrices
        auto t2 = Clock::now();
        Eigen::MatrixXi overlap = mask.array() * (*rfMask).array();
        double sum = overlap.sum();
        double dense_time = std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - t2).count();
        std::cout << "sparse time: " << sparse_time << ", dense time: " << dense_time << "\n";
        // confirm that same amount of overlap between masks is calculated
        std::cout << "sparse: " << sparse_sum << ", dense: " << sum << "\n";
        return sum;
    }

    // Given a (sparse representation) of a receptive field, check for amount of overlap and return strength of effect
    // it will have on the corresponding cell.
    double check(Eigen::SparseMatrix<int> *rfMask_sparse){
        Eigen::SparseMatrix<int> sparse_overlap = mask_sparse.cwiseProduct(*rfMask_sparse);
        double sparse_sum = sparse_overlap.nonZeros();
        return sparse_sum * amp;  // modified by intensity of stimulus
    }

    Eigen::MatrixXd getRecTable() {
        // Position, amplitude and orientation of stimulus at each time-step for this stimulus.
        // Use these to draw stimulus movies to be used by the decoder.
        int recLen = xPosRec.size();
        Eigen::MatrixXd all_recs = Eigen::MatrixXd::Ones(recLen, 5);

        for(std::size_t i = 0; i < recLen; ++i){
            all_recs(i, 0) = xPosRec[i];
            all_recs(i, 1) = yPosRec[i];
            all_recs(i, 2) = ampRec[i];
            all_recs(i, 3) = orientRec[i];
        }
        return all_recs;
    }

    void saveParams(){
        // print stim parameters to file (shape dims, orientation, angle, etc)
        // or just do this in the folder name?
    }
};


#endif //RETINA_STIM_H
