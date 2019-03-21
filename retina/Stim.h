#ifndef RETINA_STIM_H
#define RETINA_STIM_H

#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

#include "utils.h"
#include "eigen_types.h"


class Stim {
private:
    std::tuple<int, int> dims;             // dimensions of network model this cell belongs to
    Eigen::MatrixXd net_xgrid;             // X range grid used for generation of masks
    Eigen::MatrixXd net_ygrid;             // Y range grid used for generation of masks
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
    Eigen::MatrixXb mask;                  // mask defining this stimulus
    std::vector<Eigen::MatrixXb> maskRec;  // stored masks from each timestep

public:
    Stim(const std::tuple<int, int> &net_dims, const Eigen::MatrixXd &xgrid, const Eigen::MatrixXd &ygrid, const int net_dt,
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

    Eigen::MatrixXb getMask(){
        return mask;
    }

    void drawMask(){
        if (type == "bar") {
            mask = rectMask(net_xgrid, net_ygrid, pos, orient, width, length);
        } else if (type == "circle") {
            mask = circleMask(net_xgrid, net_ygrid, pos, radius);
        }
    }

    void move(){
        pos[0] += vel/dt * cos(theta_rad);
        pos[0] += vel/dt * sin(theta_rad);
        drawMask();
        maskRec.push_back(mask);
    }

    double check(Eigen::MatrixXb rfMask){
        Eigen::MatrixXb overlap = mask.array() * rfMask.array();
        return (overlap.array() == true).count();
    }
};


#endif //RETINA_STIM_H
