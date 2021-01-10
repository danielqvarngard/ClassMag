#include <iostream>
#include <fstream>

#define LATTICEFUNCTION geometry::has0

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
    auto n_overrelax = 10;
    auto n_measure = 10000;
    auto n_skip = 2;
    auto n_resamples = 100;
    unsigned int L = 4;
    const auto systemSize = std::array<unsigned int, 3>({L,L,L});
    std::array<geometry::Euclidean<3>,3> bravais;
    for (unsigned int ii = 0; ii < 3; ++ii){
        bravais[ii].fill(0.0);
        bravais[ii][ii] = 3.0;
    }
    auto lattice = geometry::Lattice<3>(bravais,systemSize);
    lattice.decorate_(geometry::icosahedralCluster());
    lattice.append_({0.0, 0.0, 0.0});
    const auto interaction = models::nearestNeighbor(-1.0,lattice,0.385);

    std::string dir = "../out/";
    std::string filename = dir + "testMC_NN_single_";
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
        1.5000, 1.2500, 1.0000, 0.8000, 0.7000, 0.6000, 0.5800, 0.5600, 0.5400, 0.5200, 0.5000, 
        0.4900, 0.4800, 0.4600, 0.4400, 0.4200, 0.4000, 0.3000, 0.2000, 0.1800, 0.1600, 0.1400,
        0.1200, 0.1100, 0.1000, 0.0900, 0.0800, 0.0600, 0.0400, 0.0200
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
        }

        auto meanEstimate = [](const std::vector<double> &x){return montecarlo::mean(x);};
        auto varianceEstimate = [](const std::vector<double> &x){return montecarlo::variance(x);};
        std::vector<double> estimates(2*n_orderParameters);

        fp << T << " ";
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