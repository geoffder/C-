#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

using namespace Eigen;


int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}

class NetworkModel {

};

class Cell {
public:
    NetworkModel model;                   // the network model this cell belongs to
    std::tuple<int, int> pos;             // centre coordinates (constant)
    int diam;                             // soma diameter
    Eigen::MatrixXi somaMask;             // mask defining cell body
    int rf;                               // receptive field radius
    Eigen::MatrixXi rfMask;               // mask defining receptive field
    float Vm = 0;                         // "membrane" state
    float dtau;                           // decay tau
    std::vector<float> rec;               // activity recording
    std::vector<std::vector<float>> recs; // collection of recordings (each trial)
};