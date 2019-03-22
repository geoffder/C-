#ifndef RETINA_STIM_H
#define RETINA_STIM_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <tuple>
#include <vector>

#include "utils.h"
#include "eigen_types.h"

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

class Stim {
private:
    int dims[2];             // dimensions of network model this cell belongs to
    Eigen::MatrixXd * net_xgrid;             // pointer to network X range grid used for generation of masks
    Eigen::MatrixXd * net_ygrid;             // pointer to network Y range grid used for generation of masks
    double dt;                              // timestep of network model
    double pos[2];                         // centre coordinates
    int tOn;
    int tOff;
    double vel;
    double theta, theta_rad;
    double orient, orient_rad;
    double amp;
    double dAmp;
    std::string type;
    double radius;                         // receptive field radius
    double length;
    double width;
    Eigen::MatrixXi mask;                  // mask defining this stimulus
    Eigen::SparseMatrix<int> mask_sparse;
    std::vector<Eigen::MatrixXi> maskRec;  // stored masks from each timestep

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

    void move(){
        pos[0] += vel/dt * cos(theta_rad);
        pos[0] += vel/dt * sin(theta_rad);
        drawMask();
        mask_sparse = mask.sparseView(1);
        maskRec.push_back(mask);
    }

    double checkOld(Eigen::MatrixXi *rfMask){
        // time comparison of mask multiplication using sparse and dense matrices
        // sparse matrices (convert first from dense)
        Eigen::SparseMatrix<int> rfMaskSparse = (*rfMask).sparseView(1);
        Eigen::SparseMatrix<int> maskSparse = mask.sparseView(1);
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

    double check(Eigen::SparseMatrix<int> *rfMask_sparse){
        Eigen::SparseMatrix<int> sparse_overlap = mask_sparse.cwiseProduct(*rfMask_sparse);
        double sparse_sum = sparse_overlap.nonZeros();
        return sparse_sum;
    }
};


#endif //RETINA_STIM_H
