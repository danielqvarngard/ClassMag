#include "Geometry/include/lattice.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;
int main(int argc, char *argv[]){

    auto n_measurements = 1000;
    auto n_thermalizations = 1000;

    std::array<geometry::Euclidean<3>,3> bravais;
    for (unsigned int ii = 0; ii < 3; ++ii){
        bravais[ii].fill(0.0);
        bravais[ii][ii] = 1.0;
    }

    const std::array<unsigned int,3> systemSize = {4, 4, 4};

    geometry::Lattice<3> lattice(bravais,systemSize);
    lattice.squareDistance_(0,1);

    auto nn_coupling = [lattice](int site1, int site2){
        if (lattice.squareDistance_(site1,site2) < 1.01 && site1 != site2)
            return 1.0;
        return 0.0;
    };
    auto seed = 0;
    montecarlo::VectorModelManager<3> mcProc(lattice.n_sites_(), nn_coupling, seed);
    mcProc.printM_();
    std::vector<double> temperature = {3.0, 2.5, 2.0, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8,
    0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06};
    std::vector<double> energy(temperature.size());
    std::vector<double> heatCapacity(temperature.size());
    auto measureindex = 0;
    auto n_sites = static_cast<double>(lattice.n_sites_());
    for (auto T : temperature){
        auto beta = 1.0/T;
        mcProc.beta_ = 1.0/T;
        std::cout << "Temperature: " << T << "\n";
        mcProc.update_(n_thermalizations);
        std::vector<double> microstateEnergy(n_measurements);
        for (unsigned int ii = 0; ii < n_measurements; ++ii){
            mcProc.update_();
            mcProc.overRelax_();
            microstateEnergy[ii] = mcProc.energy_();
        }
        auto energyEstimator = [n_sites](const std::vector<double> &x){
            return montecarlo::mean(x)/n_sites;
        };
        auto heatCapacityEstimator = [n_sites,T](const std::vector<double> &x){
            return montecarlo::variance(x)/(T*T*n_sites);
        };
        
        auto energyEstimate = montecarlo::bootstrap(microstateEnergy,energyEstimator,1000);
        auto heatCapacityEstimate = montecarlo::bootstrap(microstateEnergy,heatCapacityEstimator,1000);

        energy[measureindex] = energyEstimate.first;
        heatCapacity[measureindex] = heatCapacityEstimate.first;
        std::cout << T << " ";
        std::cout << energy[measureindex] << " ";
        std::cout << heatCapacity[measureindex] << "\n";
        ++measureindex;
        
    }
}