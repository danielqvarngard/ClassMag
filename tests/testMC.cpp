#include <iostream>

#define LATTICEFUNCTION geometry::has100

#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Interactions/include/orderParameter.hpp"
#include "Interactions/include/nearestNeighbor.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;
int main(int argc, char *argv[]){

    auto n_thermalize = 10000;
    auto n_overrelax = 10;
    auto n_measure = 1000;
    auto n_skip = 10;
    auto n_resamples = 100;

    auto const systemSize = std::array<unsigned int, 3>({4,4,4});
    auto lattice = LATTICEFUNCTION(systemSize);
    
    const auto interaction = models::nearestNeighbor(-1.0,lattice,0.385);

    auto seed = 0;
    auto mc = montecarlo::VectorModelManager<3>(
        (lattice.n_sites_()),
        interaction,
        seed);
    
    const auto mag = models::magnetization<3>();
    mc.addOrderParameter_(mag);
    if (lattice.n_decorations_() == 12){
        const auto cluster = models::clusterOrder<3,3>(lattice);
        mc.addOrderParameter_(cluster);
    }
    else if (lattice.n_decorations_() == 13){
        const auto clusterPair = models::has100Params(systemSize);
        mc.addOrderParameter_(clusterPair.first);
        mc.addOrderParameter_(clusterPair.second);
    }

    
    
    auto temperatures = {
        3.0, 2.5, 2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6,
        0.58, 0.56, 0.54, 0.52, 0.5, 0.48, 0.46, 0.44, 0.42, 0.4,
        0.3, 0.2, 0.1, 0.08, 0.06, 0.04};

    auto heatCapMeasures = std::vector<double>();
    auto magMeasures = std::vector<double>();
    auto chiMeasures = std::vector<double>();

    for (auto T : temperatures){
        mc.beta_ = 1.0/T;
        mc.update_(n_thermalize);
        std::vector<double> energyData(n_measure);
        std::vector<double> magData(n_measure);
        std::vector<double> clusterData(n_measure);
        for (unsigned int ii = 0; ii < n_measure; ++ii){
            mc.overRelax_(n_overrelax);
            mc.update_(n_skip);
            energyData[ii] = mc.energy_();
            auto mag = mc.measure_();
            magData[ii] = mag[0];
            clusterData[ii] = mag[1];
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

        auto clusterResult = montecarlo::bootstrap(clusterData, meanEstimate, n_resamples);

        heatCapMeasures.push_back(heatCapResult.first);
        magMeasures.push_back(magResult.first);
        chiMeasures.push_back(chiResult.first);

        std::cout << T << " ";
        std::cout << montecarlo::mean(energyData)/static_cast<double>(lattice.n_sites_()) << " ";
        std::cout << heatCapResult.first << " ";
        std::cout << magResult.first << " ";
        std::cout << clusterResult.first << " ";
        std::cout << chiResult.first << "\n";
    }
    
}