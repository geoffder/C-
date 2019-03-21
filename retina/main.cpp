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
    std::cout << "Hello, World!" << std::endl;

    int dims[2] = {20, 20};
    double dt = 1;
    double pos[2] = {10, 10};
    double diam = 5;
    double rf = 7;
    double dtau = 1;


    auto [xgrid, ygrid] = gridMats(dims[0], dims[1]);

    Cell a_cell(dims, xgrid, ygrid, dt, pos, diam, rf, dtau);

    std::cout << a_cell.getSoma() << "\n\n" << a_cell.getRF() << "\n\n";

    a_cell.setVm(10);
    std::cout << "voltage: " << a_cell.getVm() << "\n";

    for(int i = 0; i < 10; ++i){
      a_cell.decay();
    }

    std::vector<double> recording = a_cell.getRec();

    for(int i = 0; i < 10; ++i){
        std::cout << recording[i] << " ";
    }
    std::cout << "\n";

    return 0;
}