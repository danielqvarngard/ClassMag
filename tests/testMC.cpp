#include <iostream>
#include <fstream>

#define LATTICEFUNCTION geometry::has100

#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Interactions/include/orderParameter.hpp"
#include "Interactions/include/nearestNeighbor.hpp"
#include "Interactions/include/rkky.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;
int main(int argc, char *argv[]){

    auto n_thermalize = 10000;
    auto n_overrelax = 5;
    auto n_measure = 10000;
    auto n_skip = 2;
    auto n_resamples = 100;
    unsigned int L = 6;
    auto const systemSize = std::array<unsigned int, 3>({L,L,L});
    auto lattice = LATTICEFUNCTION(systemSize);
    int kf_label = 18;
    auto kf = static_cast<double>(kf_label);
    const auto interaction = models::rkkyInteraction(kf,lattice,4.0);

    std::string dir = "../out/";
    std::string filename = dir + "testMC_RKKY_";
    filename += std::to_string(kf_label);
    filename += "_";
    auto seed = 0;
    auto mc = montecarlo::VectorModelManager<3>(
        (lattice.n_sites_()),
        interaction,
        seed);
    
    const auto mag = models::magnetization<3>();
    mc.addOrderParameter_(mag);
    auto n_orderParameters = 1;
    
    if (lattice.n_decorations_() == 12){
        filename += "has0_";
        const auto cluster = models::clusterOrder<3,3>(lattice);
        mc.addOrderParameter_(cluster);
        ++n_orderParameters;
    }
    else if (lattice.n_decorations_() == 13){
        filename += "has100_";
        const auto clusterPair = models::has100Params(systemSize);
        mc.addOrderParameter_(clusterPair.first);
        mc.addOrderParameter_(clusterPair.second);
        n_orderParameters = n_orderParameters + 2;
    }
    filename += "L_";
    filename += std::to_string(L);
    filename += ".mcout";

    std::vector<std::vector<double>> opc(n_orderParameters);
    for (unsigned int jj = 0; jj < n_orderParameters; ++jj)
        opc[jj].resize(n_measure);
    
    auto temperatures = {
        0.0005, 0.0004, 0.00035, 0.000325, 0.0003, 0.000275, 0.00025, 0.000225, 0.0002, 0.000175,
        0.00015, 0.000125, 0.0001
    };

    auto fp = std::ofstream(filename);
    auto zeroPatience = 1;
    for (auto T : temperatures){
        mc.beta_ = 1.0/T;
        std::vector<double> energy(n_measure);
        mc.update_(n_thermalize,n_overrelax);
        for (unsigned int ii = 0; ii < n_measure; ++ii){
            mc.update_(n_skip,n_overrelax);
            energy[ii] = mc.energy_();
            auto mag = mc.measure_();
            for (unsigned int jj = 0; jj < n_orderParameters; ++jj){
                opc[jj][ii] = mag[jj];
            }
        }

        auto meanEstimate = [](const std::vector<double> &x){return montecarlo::mean(x);};
        auto varianceEstimate = [](const std::vector<double> &x){return montecarlo::variance(x);};
        std::vector<double> estimates(2*n_orderParameters);

        fp << T << " ";
        for (unsigned int jj = 0; jj < n_orderParameters; ++jj){
            auto meanval = montecarlo::bootstrap(opc[jj],meanEstimate,n_resamples);
            auto n_sites = static_cast<double>(lattice.n_sites_());
            fp << meanval.first/n_sites << " ";
            auto varval = montecarlo::bootstrap(opc[jj],varianceEstimate,n_resamples);
            fp << varval.first/(T*n_sites) << " ";
        }
        {
        auto meanval = montecarlo::bootstrap(energy,meanEstimate,n_resamples);
        auto n_sites = static_cast<double>(lattice.n_sites_());
        fp << meanval.first/n_sites << " ";
        fp << meanval.second/n_sites << " ";
        auto varval = montecarlo::bootstrap(energy,varianceEstimate,n_resamples);
        fp << varval.first/(T*T*n_sites) << " ";
        fp << varval.second/(T*T*n_sites);
        }
        fp << "\n";
        std::cout << zeroPatience << "/" << temperatures.size();
        std::cout << " temperatures completed\n";
        ++zeroPatience;

    }

    fp.close();
    
}