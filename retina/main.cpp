#include <windows.h>
#include <iostream>
#include <tuple>
#include <vector>
//#include <filesystem>

#include <Eigen/Dense>

//#include <boost/filesystem.hpp>
//#define BOOST_NO_CXX11_SCOPED_ENUMS
//#include <boost/filesystem.hpp>
//#undef BOOST_NO_CXX11_SCOPED_ENUMS

#include "type_defs.h"
#include "utils.h"

#include "NetworkModel.h"
#include "Cell.h"
#include "Stim.h"


using namespace Eigen;
//namespace fs = boost::filesystem;
//namespace fs = std::filesystem;

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

    int dims[2] = {20, 20};
    double dt = 1;
    double pos[2] = {10, 10};
    double diam = 5;
    double rf = 7;
    double dtau = 1;


    auto [xgrid, ygrid] = gridMats(dims[0], dims[1]);

    Cell a_cell(dims, xgrid, ygrid, dt, pos, diam, rf, dtau);

    std::cout << a_cell.getSoma() << "\n\n" << a_cell.getRF() << "\n\n";
    MatrixXi soma = a_cell.getSoma();

    a_cell.setVm(10);
    std::cout << "voltage: " << a_cell.getVm() << "\n\n";

    for(int i = 0; i < 10; ++i){
      a_cell.decay();
    }

    std::vector<double> recording = a_cell.getRec();

    for(int i = 0; i < 10; ++i){
        std::cout << recording[i] << " ";
    }
    std::cout << "\n";
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
        net.newStim(start_pos, 0, 500, double(1), dir, -dir, 1, 0, "bar", 0, 50, 100);
        // net.newStim(start_pos, 0, 500, double(1), dir, -dir, 1, 0, "circle", 50);
        net.run();
        net.clearStims();
    }
    std::cout << std::endl;

    std::cout << "making network matrix...\n\n";
    MatrixXd cellMat = net.cellMatrix();

    // make folder (Windows), fails if folder already exists without issue
    CreateDirectory((baseFolder + dataPrefix).c_str(), nullptr);

    MatrixXdToCSV(baseFolder + dataPrefix + "/cellMat.csv", cellMat);
    MatrixXdToCSV(baseFolder + dataPrefix + "/cellCoords.csv", net.getCellXYs());
    MatrixXdToCSV(baseFolder + dataPrefix + "/cellRecs.csv", net.getRecTable());

    // wait for ENTER to close terminal
    std::string wait;
    std::cout << "Type anything and hit ENTER to terminate..." << std::endl;
    std::cin >> wait;

    return 0;
}