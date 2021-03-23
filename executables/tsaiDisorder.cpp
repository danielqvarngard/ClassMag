#include "Geometry/include/hasx.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/rkky.hpp"
#include "Base/include/hasxStagMag.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/mcProfile.hpp"
#include "MonteCarlo/include/postProcessing.hpp"
#include <array>
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <mpi.h>

using namespace classmag;

int main(int argc, char** argv){
    MPI_Init(NULL,NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    const int percentage = 63u;
    auto occupancy = static_cast<double>(percentage)/100.0;
    const auto L = 4u;
    const std::array<unsigned int, 3> size({L,L,L});
    const auto lattice = geometry::hasx(size, occupancy, world_rank);
    auto partitions = lattice.partitions_();
    auto n_shell = static_cast<double>(partitions[2]);
    auto n_center = static_cast<double>(lattice.n_sites_()) - n_shell;

    std::string dir = "../out/";
    std::string filename = dir;
    filename += "has" + std::to_string(percentage);
    filename += "_L_" + std::to_string(L);
    filename += "_proc_" + std::to_string(world_rank);
    filename += ".mcout";

    montecarlo::VectorModel_Profile mcp;
    auto n_resamples = 10;
    mcp.measurement_ = 10000;
    mcp.thermalization_ = 10000;
    mcp.skips_ = 1;
    mcp.n_sites_ = lattice.n_sites_();
    mcp.overrelax_ = 0;
    mcp.seed_ = 138;
    mcp.getPartitions_(lattice);
    auto n_sublattices = mcp.partitions_.size() - 1;

    const auto interaction = base::rkkyInteraction(19.8158,lattice, 2.0*L);
    auto mc = montecarlo::VectorModelManager<3>(
        mcp,
        interaction);

    auto opPair = base::chasxMag_split(lattice);

    mc.addOrderParameter_(opPair.first);
    mc.addOrderParameter_(opPair.second);

    auto temperatures = {
        0.001, 0.0008, 0.0006, 0.0005, 0.00045, 0.0004, 0.00035, 0.0003, 0.00025
    };

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

    return 0;
}