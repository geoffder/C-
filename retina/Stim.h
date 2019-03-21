#ifndef RETINA_STIM_H
#define RETINA_STIM_H

//#include <iostream>
//#include <Eigen/Dense>
//#include <tuple>
//#include <vector>
//
//#include "utils.h"
//#include "eigen_types.h"
#include "everything.h"
//using namespace utils;
using namespace Eigen;
//typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

//std::tuple<MatrixXd, MatrixXd> rotateGrids(double origin[2], MatrixXd &xgrid, MatrixXd &ygrid, double degrees);
//MatrixXb circleMask(MatrixXd &xgrid, MatrixXd &ygrid, double origin[2], double radius);
//MatrixXb rectMask(MatrixXd &xgrid, MatrixXd &ygrid, double origin[2], double orient, double width, double height);

class Stim {
private:
    std::tuple<int, int> dims;            // dimensions of network model this cell belongs to
    MatrixXd net_xgrid;                   // X range grid used for generation of masks
    MatrixXd net_ygrid;                   // Y range grid used for generation of masks
    float dt;                             // timestep of network model
    double pos[2];                        // centre coordinates
    int tOn;
    int tOff;
    double vel;
    double theta, theta_rad;
    double orient, orient_rad;
    double amp;
    double dAmp;
    std::string type;
    double radius;                        // receptive field radius
    double length;
    double width;
    MatrixXb mask;                        // mask defining this stimulus
    std::vector<MatrixXb> maskRec;        // stored masks from each timestep

public:
    Stim(const std::tuple<int, int> &net_dims, const MatrixXd &xgrid, const MatrixXd &ygrid, const int net_dt,
            const int start_pos[2], const int time_on, const int time_off, const double velocity,
            const double direction, const double orientation, const double amplitude, const double change){
        // network properties
        dims = net_dims;
        net_xgrid = xgrid;
        net_ygrid = ygrid;
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

    MatrixXb getMask(){
        return mask;
    }

    void drawMask(){
        if (type == "bar") {
            mask = utils::rectMask(net_xgrid, net_ygrid, pos, orient, width, length);
        } else if (type == "circle") {
            mask = utils::circleMask(net_xgrid, net_ygrid, pos, radius);
        }
    }

    void move(){
        pos[0] += vel/dt * cos(theta_rad);
        pos[0] += vel/dt * sin(theta_rad);
        drawMask();
        maskRec.push_back(mask);
    }

    float check(MatrixXb rfMask){
        MatrixXb overlap = mask.array() * rfMask.array();
        return (overlap.array() == true).count();
    }
};


#endif //RETINA_STIM_H
