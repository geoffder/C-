#include <iostream>
#include <Eigen/Dense>
#include <tuple>
#include <vector>

using namespace Eigen;
typedef Matrix<bool, Dynamic, Dynamic> MatrixXb;

class NetworkModel {

};

class Cell {
private:
    std::tuple<int, int> dims;            // dimensions of network model this cell belongs to
    int dt;                               // timestep of network model
    std::tuple<int, int> pos;             // centre coordinates (constant)
    int diam;                             // soma diameter
    Eigen::MatrixXi somaMask;             // mask defining cell body
    int rf;                               // receptive field radius
    Eigen::MatrixXi rfMask;               // mask defining receptive field
    float Vm;                             // "membrane" state
    float dtau;                           // decay tau
    std::vector<float> rec;               // activity recording
    std::vector<std::vector<float>> recs; // collection of recordings (each trial)

public:
    Cell(){
        Vm = 0;
    }
    void setParams(const std::tuple<int, int> &net_dims, int net_dt, const std::tuple<int, int> &cell_pos,
                   int cell_diam, int rf_rad, float cell_dtau) {
        dims = net_dims;
        dt = net_dt;
        pos = cell_pos;
        diam = cell_diam;
        rf = rf_rad;
        dtau = cell_dtau;
    }

    float getVm(){
        return Vm;
    }

    std::vector<float> getRec(){
        return rec;
    }

    void setVm(float newVm){
        Vm = newVm;
    }

    void storeRec(){
        recs.push_back(rec);
        rec.clear();
    }

    void excite(float strength){
        Vm += strength;
    }

    void decay(){
        float delta = Vm * (1 - exp(-dt/dtau));
        Vm = std::fmax(float (0), Vm - delta);
        rec.push_back(Vm);
    }
};

int main() {
    std::cout << "Hello, World!" << std::endl;

    std::tuple<int, int> dims = std::make_tuple(600, 600);
    int dt = 1;
    std::tuple<int, int> pos = std::make_tuple(50, 50);
    int diam = 10;
    int rf = 50;
    float dtau = 1;

    Cell a;
    a.setParams(dims, dt, pos, diam, rf, dtau);
    a.setVm(10);
    std::cout << "voltage: " << a.getVm() << "\n";

    for(int i = 0; i < 10; ++i){
      a.decay();
    }

    std::vector<float> recording = a.getRec();

    for(int i = 0; i < 10; ++i){
        std::cout << recording[i] << " ";
    }
    std::cout << "\n";

    return 0;
}