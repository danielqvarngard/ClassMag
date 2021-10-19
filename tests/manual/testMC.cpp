#include <iostream>
#include <fstream>

#define LATTICEFUNCTION geometry::has0

#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Base/include/orderParameter.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/rkky.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;
int main(int argc, char *argv[]){
    montecarlo::VectorModel_Profile mcp;
    auto n_thermalize = 10000;
    auto n_overrelax = 1;
    auto n_measure = 10000;
    auto n_skip = 1;
    auto n_resamples = 100;
    unsigned int L = 5;
    const auto systemSize = std::array<unsigned int, 3>({L,L,L});

    #if 0
    /* bcc lattice for testing: */
    auto sublattice1 = geometry::cubicLattice<3>(systemSize);
    auto lattice = geometry::Lattice<3>(sublattice1);
    auto sublattice2 = geometry::cubicLattice<3>(systemSize);
    auto e = geometry::Euclidean<3>({0.5, 0.5, 0.5});
    sublattice2.decorate_({e});

    lattice.append_(sublattice2);
    const auto interaction = base::nearestNeighbor(-1.0,lattice,0.9);
    #endif
    
    
    auto lattice = geometry::chas100(systemSize);
    const auto interaction = base::nearestNeighbor(-1.0,lattice,0.385);
   
    mcp.measurement_ = n_measure;
    mcp.thermalization_ = n_thermalize;
    mcp.skips_ = n_skip;
    mcp.n_sites_ = lattice.n_sites_();
    mcp.overrelax_ = n_overrelax;
    mcp.seed_ = 138;
    mcp.getPartitions_(lattice);
    auto n_sublattices = mcp.partitions_.size() - 1;

    std::cout << mcp.n_sites_ << "\n";

    std::string dir = "../out/";
    std::string filename = dir + "testMC_testmag_";
    auto mc = montecarlo::VectorModelManager<3>(
        mcp,
        interaction);
    
    filename += "L_bugtest";
    filename += std::to_string(L);
    filename += ".mcout";
    
    #if 1
    auto temperatures = {
        1.2500, 1.0000, 0.8000, 0.7000, 0.6000, 0.5800, 0.5600, 0.5400, 0.5200, 0.5000, 
        0.4900, 0.4800, 0.4600, 0.4400, 0.4200, 0.4000, 0.3000, 0.2000, 0.1800, 0.1600, 0.1400,
        0.1200, 0.1100, 0.1000, 0.0900, 0.0800, 0.0600, 0.0400
    };
    #endif

    #if 0
    auto temperatures = {1.5, 1.25, 1.0};
    #endif

    auto fp = std::ofstream(filename);
    auto zeroPatience = 1;
    for (auto T : temperatures){
        mc.beta_ = 1.0/T;
        mc.update_(n_thermalize,n_overrelax);
        std::vector<base::MagnetizationData<3,2>> magdat(n_measure);
        for (unsigned int ii = 0; ii < n_measure; ++ii){
            mc.update_(n_skip,n_overrelax);
            magdat[ii] = mc.magnetization_();
        }
        fp << T << " ";
        auto sus = montecarlo::trChiBootstrap(magdat, 1.0/T, n_resamples);
        for (auto ii = 0u; ii < n_sublattices; ++ii){
            fp << sus[ii].first;
            fp << " ";
            fp << sus[ii].second;
            fp << " ";
        }
        fp << "\n";
        std::cout << zeroPatience << "/" << temperatures.size();
        std::cout << " temperatures completed\n";
        ++zeroPatience;

    }

    fp.close();
    
}