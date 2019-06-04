#include <windows.h>
#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
//#include <eigen3-hdf5.hpp>

#include "type_defs.h"
#include "utils.h"

#include "NetworkModel.h"
#include "Cell.h"
#include "Stim.h"


using namespace Eigen;

int main() {
    // Eigen::setNbThreads(2);
    std::cout << "threads used by Eigen: " << Eigen::nbThreads( ) << std::endl;
    std::string baseFolder = "D://retina-sim-data/";
    std::string netFolder;
    std::cout << "Use base folder " << baseFolder << "? [Y/N]" << std::endl;
    std::string answer;
    std::cin >> answer;
    if (answer.find('N') != std::string::npos || answer.find('n') != std::string::npos){
        std::cout << "Enter new base folder path:" << std::endl;
        std::cin >> baseFolder;
    }

    // std::string stim_type[4] = {"light_bar", "light_circle", "dark_bar", "dark_circle"};
    // std::string stim_type[4] = {"thin_light_bar", "thin_dark_bar", "thick_light_bar", "thick_dark_bar"};
    std::string stim_type[6] = {"thin_light_bar", "med_light_bar", "thick_light_bar", "thin_dark_bar", "med_dark_bar", "thick_dark_bar"};
    // std::string stim_type[4] = {"thin_light_ellipse", "thin_dark_ellipse", "thick_light_ellipse", "thick_dark_ellipse"};
    double directions[8] = {0, 45, 90, 135, 180, 225, 270, 315};
    int net_dims[2] = {700, 700};
    NetworkModel net(net_dims, int(200), int(3000), double(5));

    for(int i = 0; i < 20; ++i) {
        std::cout << "Constructing net" << i << "..." << std::endl;
        net.populate(int(20), double(10));
        std::cout << "Number of cells: " << net.getCells().size() << std::endl;


        netFolder = baseFolder+"net"+std::to_string(i)+"/";
        CreateDirectory((netFolder).c_str(), nullptr);
        for (auto &stm : stim_type) {
            std::cout << "Running " << stm << "... \ndirs:" << std::endl;
            for (auto &dir : directions) {
                std::cout << dir << " ";
                auto[cx, cy] = net.getOrigin();
                double start_pos[2] = {cx - cx * cos(deg2rad(dir))*1.3, cy - cy * sin(deg2rad(dir))*1.3};

                if (stm == "thin_light_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "bar", 0, 30, 600);
                } else if (stm == "med_light_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "bar", 0, 100, 600);
                } else if (stm == "thick_light_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "bar", 0, 400, 600);
                } else if (stm == "thin_light_ellipse") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "ellipse", 0, 30, 600);
                } else if (stm == "thick_light_ellipse") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "ellipse", 0, 200, 600);
                } else if (stm == "light_circle") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, 1, 0, "circle", 100);
                } else if (stm == "thin_dark_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "bar", 0, 30, 600);
                } else if (stm == "med_dark_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "bar", 0, 100, 600);
                } else if (stm == "thick_dark_bar") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "bar", 0, 400, 600);
                } else if (stm == "thin_dark_ellipse") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "ellipse", 0, 30, 600);
                } else if (stm == "thick_dark_ellipse") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "ellipse", 0, 200, 600);
                } else if (stm == "dark_circle") {
                    net.newStim(start_pos, 0, 3000, double(.6), dir, -dir, -1, 0, "circle", 100);
                }
                net.run(netFolder, stm + std::to_string(std::lround(dir)));
                net.clearStims();
            }
        }
        std::cout << "saving network information...\n\n";
        net.netToFile(netFolder);
        net.clearCells();
    }

    // wait for ENTER to close terminal
    std::string wait;
    std::cout << "Type anything and hit ENTER to terminate..." << std::endl;
    std::cin >> wait;

    return 0;
}