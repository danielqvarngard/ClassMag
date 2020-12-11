#include <iostream>

#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Interactions/include/orderParameter.hpp"
#include "Interactions/include/nearestNeighbor.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;
int main(int argc, char *argv[]){

    auto n_thermalize = 1000;
    auto n_overrelax = 100;
    auto n_measure = 1000;
    auto n_resamples = 1000;

    auto const systemSize = std::array<unsigned int, 3>({4,4,4});
    auto lattice = geometry::cubicLattice<3>(systemSize);
    const auto interaction = models::nearestNeighbor(1.0,lattice,1.001);
    const auto mag = models::magnetization<3>();
    auto seed = 0;
    auto mc = montecarlo::VectorModelManager<3>(
        (lattice.n_sites_()),
        interaction,
        seed);

    mc.addOrderParameter_(mag);

    auto temperatures = {
        3.0, 2.5, 2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.05, 1.025, 1.0, 
        0.975, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4,
        0.3, 0.2, 0.1, 0.08, 0.06, 0.04};

    auto heatCapMeasures = std::vector<double>();
    auto magMeasures = std::vector<double>();
    auto chiMeasures = std::vector<double>();

    for (auto T : temperatures){
        mc.beta_ = 1.0/T;
        mc.update_(n_thermalize);
        std::vector<double> energyData(n_measure);
        std::vector<double> magData(n_measure);
        for (unsigned int ii = 0; ii < n_measure; ++ii){
            mc.overRelax_(n_overrelax);
            mc.update_();
            energyData[ii] = mc.energy_();
            auto mag = mc.measure_();
            magData[ii] = mag[0];
        }

        auto meanEstimate = [](const std::vector<double> &x){return montecarlo::mean(x);};
        auto varianceEstimate = [](const std::vector<double> &x){return montecarlo::variance(x);};

        auto heatCapResult = montecarlo::bootstrap(energyData,varianceEstimate,n_resamples);
        heatCapResult.first /= T*T;
        heatCapResult.first /= static_cast<double>(lattice.n_sites_());

        auto magResult = montecarlo::bootstrap(magData,meanEstimate,n_resamples);
        magResult.first /= static_cast<double>(lattice.n_sites_());
        auto chiResult = montecarlo::bootstrap(magData,varianceEstimate,n_resamples);
        chiResult.first /= static_cast<double>(lattice.n_sites_());
        chiResult.first /= T;

        heatCapMeasures.push_back(heatCapResult.first);
        magMeasures.push_back(magResult.first);
        chiMeasures.push_back(chiResult.first);

        std::cout << T << " ";
        std::cout << montecarlo::mean(energyData)/static_cast<double>(lattice.n_sites_()) << " ";
        std::cout << heatCapResult.first << " ";
        std::cout << magResult.first << " ";
        std::cout << chiResult.first << "\n";
    }
    
}