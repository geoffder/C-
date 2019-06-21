#include <iostream>
#include <fstream>
#include <tuple>
#include <vector>
#include <math.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "type_defs.h"
#include "utils.h"
#include "stimuli.h"

#include "NetworkModel.h"

using namespace Eigen;

typedef void (*Stimuli) (NetworkModel &net, double cx, double cy, std::string &netFolder);

double DIRECTIONS[8] = {0, 45, 90, 135, 180, 225, 270, 315};

std::array<double, 2> getStartPos(double cx, double cy, double dir){
    return std::array<double, 2> ({cx - cx * cos(deg2rad(dir)) * 2, cy - cy * sin(deg2rad(dir)) * 2});
}

void dir_thin_light_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thin_light_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "bar", 0, 15, 300);
        net.run(netFolder, "thin_light_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_med_light_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running med_light_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "bar", 0, 50, 300);
        net.run(netFolder, "med_light_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thick_light_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thick_light_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "bar", 0, 200, 300);
        net.run(netFolder, "thick_light_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thin_dark_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thin_dark_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "bar", 0, 15, 300);
        net.run(netFolder, "thin_dark_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_med_dark_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running med_dark_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "bar", 0, 50, 300);
        net.run(netFolder, "med_dark_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thick_dark_bar(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thick_dark_bar... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "bar", 0, 200, 300);
        net.run(netFolder, "thick_dark_bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thin_light_ellipse(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thin_light_ellipse... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "ellipse", 0, 15, 300);
        net.run(netFolder, "thin_light_ellipse" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thick_light_ellipse(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thick_light_ellipse... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "ellipse", 0, 100, 300);
        net.run(netFolder, "thick_light_ellipse" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thin_dark_ellipse(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thin_dark_ellipse... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "ellipse", 0, 15, 300);
        net.run(netFolder, "thin_dark_ellipse" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_thick_dark_ellipse(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running thick_dark_ellipse... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "ellipse", 0, 100, 300);
        net.run(netFolder, "thick_dark_ellipse" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_small_light_circle(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running small_light_circle... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, 1, 0, "circle", 50);
        net.run(netFolder, "small_light_circle" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void dir_small_dark_circle(NetworkModel &net, double cx, double cy, std::string &netFolder) {
    std::cout << "Running small_dark_circle... \ndirs:" << std::endl;
    for (auto &dir : DIRECTIONS) {
        std::cout << dir << " ";
        auto start_pos = getStartPos(cx, cy, dir);
        net.newStim(start_pos, 0, 3000, double(.3), dir, -dir, -1, 0, "circle", 50);
        net.run(netFolder, "small_dark_circle" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
}

void runExperiment_1(NetworkModel &net, std::string &netFolder) {
    auto[cx, cy] = net.getOrigin();

    std::vector<Stimuli> stim_list = {
            dir_thin_light_bar,
            dir_med_light_bar,
            dir_thick_light_bar,
            dir_thin_dark_bar,
            dir_med_dark_bar,
            dir_thick_dark_bar,
            dir_thin_light_ellipse,
            dir_thick_light_ellipse,
            dir_thin_dark_ellipse,
            dir_thick_dark_ellipse,
            dir_small_light_circle,
            dir_small_dark_circle

    };

    for (auto &stimulus : stim_list){
        stimulus(net, cx, cy, netFolder);
    }
}