#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

using Eigen;


int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}

class Cell {
public:
    NetworkModel model;  // the network model this cell belongs to
    std::tuple<int, int> pos;  // centre coordinates (constant)
    int diam;  // soma diameter
    Eigen::MatrixXi somaMask;
    int rf  // receptive field radius
    Eigen::MatrixXi rfMask  //
    float Vm = 0
    float dtau  // decay tau
    std::vector<float> rec
    std::vector<std::vector<float>> recs
};