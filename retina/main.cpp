#include <windows.h>
#include <iostream>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
//#include <eigen3-hdf5.hpp>

#include "type_defs.h"
#include "utils.h"
#include "stimuli.h"

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


    std::array<int, 2> net_dims = {350, 350};  // 700, 700
    NetworkModel net(net_dims, int(100), int(3000), double(5));  // 200 margins for 700, 700

    for(int i = 0; i < 20; ++i) {
        // fill the empty network object with cells
        std::cout << "Constructing net" << i << "..." << std::endl;
        net.populate(int(10), double(5));
        std::cout << "Number of cells: " << net.getCells().size() << std::endl;

        // create network directory
        netFolder = baseFolder+"net"+std::to_string(i)+"/";
        CreateDirectory((netFolder).c_str(), nullptr);

        // run series of stimuli, saving results along the way
        runExperiment_1(net, netFolder);

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