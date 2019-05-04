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
    std::cout << "Use base folder " << baseFolder << "? [Y/N]" << std::endl;
    std::string answer;
    std::cin >> answer;
    if (answer.find('N') != std::string::npos || answer.find('n') != std::string::npos){
        std::cout << "Enter new base folder path:" << std::endl;
        std::cin >> baseFolder;
    }

    std::string dataPrefix = "run";

    int net_dims[2] = {600, 600};
    NetworkModel net(net_dims, int(100), int(500), double(1));
    net.populate(int(40), double(10));

    std::cout << "Number of cells: " << net.getCells().size() << "\n\n";

    double directions[8] = {0, 45, 90, 135, 180, 225, 270, 315};
    std::cout << "Running... \ndirs:" << std::endl;
    for (auto& dir : directions) {
        std::cout << dir << " ";
        auto [cx, cy] = net.getOrigin();
        double start_pos[2] = {cx - cx*cos(deg2rad(dir)), cy - cy*sin(deg2rad(dir))};
        net.newStim(start_pos, 0, 500, double(1), dir, -dir, 1, 0, "bar", 0, 50, 400);
        // net.newStim(start_pos, 0, 500, double(1), dir, -dir, 1, 0, "circle", 50);
        CreateDirectory((baseFolder+"net0/").c_str(), nullptr);
        net.run(baseFolder+"net0/", "bar" + std::to_string(std::lround(dir)));
        net.clearStims();
    }
    std::cout << std::endl;

    std::cout << "saving network information...\n\n";
    net.netToFile(baseFolder+"net0/");

    // wait for ENTER to close terminal
    std::string wait;
    std::cout << "Type anything and hit ENTER to terminate..." << std::endl;
    std::cin >> wait;

    return 0;
}