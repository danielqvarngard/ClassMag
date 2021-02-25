#include "Geometry/include/hasx.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/hasxStagMag.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/mcProfile.hpp"
#include "MonteCarlo/include/postProcessing.hpp"
#include <array>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace classmag;

int main(){
    auto percentage = 100u;
    auto occupancy = static_cast<double>(percentage)/100.0;
    auto L = 4u;
    const std::array<unsigned int,3> size({L,L,L});
    auto lattice = geometry::hasx(size, occupancy);
    auto partitions = lattice.partitions_();
    auto n_shell = static_cast<double>(partitions[2]);
    auto n_center = static_cast<double>(lattice.n_sites_()) - n_shell;
    
    montecarlo::VectorModel_Profile mcp;
    auto n_resamples = 100;
    mcp.measurement_ = 50000;
    mcp.thermalization_ = 50000;
    mcp.skips_ = 1;
    mcp.n_sites_ = lattice.n_sites_();
    mcp.overrelax_ = 2;
    mcp.seed_ = 138;
    mcp.getPartitions_(lattice);
    auto n_sublattices = mcp.partitions_.size() - 1;

    const auto interaction = base::nearestNeighbor(-1.0,lattice,0.385);
    auto mc = montecarlo::VectorModelManager<3>(
        mcp,
        interaction);

    auto opPair = base::chasxStagMag_split(lattice);

    mc.addOrderParameter_(opPair.first);
    mc.addOrderParameter_(opPair.second);
    
    
    std::cout << mcp.n_sites_ << "\n";

    std::string dir = "../out/";
    std::string filename = dir + "stagmag_";
    
    filename += "L_";
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
        mc.thermalize_();
        std::vector<double> shellData(mcp.measurement_);
        std::vector<double> centerData(mcp.measurement_);
        for (unsigned int ii = 0; ii < mcp.measurement_; ++ii){
            mc.update_(mcp.skips_,mcp.overrelax_);
            auto data = mc.measure_();
            shellData[ii] = data[0];
            centerData[ii] = data[1];
        }
        fp << T << " ";
        auto shellResult = montecarlo::defaultBootstrap(shellData, n_resamples);
        auto centerResult = montecarlo::defaultBootstrap(centerData, n_resamples);
        fp << shellResult.meanEstimate/n_shell << " " << shellResult.meanDeviation/n_shell << " ";
        fp << shellResult.varianceEstimate/(T * n_shell) << " ";
        fp << shellResult.varianceDeviation/(T * n_shell) << " ";
        fp << centerResult.meanEstimate/n_center << " " << centerResult.meanDeviation/(n_center) << " ";
        fp << centerResult.varianceEstimate/(T * n_center) << " ";
        fp << centerResult.varianceDeviation/(T * n_center) << " ";
        fp << "\n";
        std::cout << zeroPatience << "/" << temperatures.size();
        std::cout << " temperatures completed\n";
        ++zeroPatience;

    }

    fp.close();
};