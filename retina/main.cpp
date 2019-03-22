#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

#include "NetworkModel.h"
#include "Cell.h"
#include "Stim.h"
#include "eigen_types.h"
#include "NetworkModel.h"
#include "utils.h"

using namespace Eigen;

int main() {
//    Eigen::setNbThreads(2);
    std::cout << "threads used by Eigen: " << Eigen::nbThreads( ) << std::endl;
    std::string dataFolder = "D://work/";

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

    double start_pos[2] = {300, 300};
    // defaults don't work the same in C++ as python, must to arguments left to right, no skipping.
    // so if I do a bar, I still have to provide something for radius. For circles I can end with radius (ignore WxL).
    net.newStim(start_pos, 0, 500, double(1), double(0), double(0), double(1), double(0), "bar", 0, 50, 100);
//    net.newStim(start_pos, 0, 500, double(1), double(0), double(0), double(1), double(0), "circle", 50);
    auto stims = net.getStims();
    stims[0].drawMask();
    MatrixXi stimMask = stims[0].getMask();

//    auto [xgrid_full, ygrid_full] = gridMats(600, 600);
//    MatrixXi bigRF = circleMask(xgrid_full, ygrid_full, start_pos, 50);
//    double strength = stims[0].check(bigRF);
//    std::cout << "check result: " << strength;

    net.run();

    std::cout << "making network matrix...\n\n";
    MatrixXd cellMat = net.cellMatrix();
    std::cout << "writing files...\n\n";
    MatrixXiToCSV(dataFolder + "mask.csv", soma);
    MatrixXdToCSV(dataFolder + "cellMat.csv", cellMat);
    MatrixXdToCSV(dataFolder + "cellCoords.csv", net.getCellXYs());
    MatrixXiToCSV(dataFolder + "stimMask.csv", stimMask);
    MatrixXdToCSV(dataFolder + "cellRecs.csv", net.getRecTable());
    return 0;
}