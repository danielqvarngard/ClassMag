#include "Base/include/nearestNeighbor.hpp"
#include "FileIO/include/StreamManager.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/extendedEnsembles.hpp"
using namespace classmag;
int main(int argc, char* argv[]){
    auto str = fileio::datestamp();
    auto filename = str + ".pttest";
    auto os = fileio::OStreamManager(filename);
    auto sweeps = 1u;
    const auto temperatures = {2.0, 1.5, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.08};
    auto betas = std::vector<double>();
    for (auto T : temperatures)
        betas.push_back(1.0/T);   
    os << betas;
    auto L = 4u;
    const auto size = std::array<unsigned int,3>({L,L,L});
    auto sublattice = geometry::cubicLattice<3>(size);
    auto lattice = geometry::Lattice(sublattice);
    auto interaction = base::nearestNeighbor(-1.0,lattice,1.01);
    auto pt = montecarlo::ParallelTemperer(betas);
    auto mcp = montecarlo::VectorModel_Profile();
    mcp.measurement_ = 10000;
    mcp.thermalization_ = 10000;
    mcp.skips_ = 1;
    mcp.overrelax_ = 2;
    mcp.n_sites_ = lattice.n_sites_();
    auto mc = montecarlo::VectorModelManager<3>(mcp,interaction);
    auto mcs = std::vector<montecarlo::VectorModelManager<3>>(betas.size(),mc);

    for (auto ii = 0u; ii < mcs.size(); ++ii){
        mcs[ii].seed_(ii);
        mcs[ii].beta_ = betas[ii];
    }

    for (auto thermsteps = 0u; thermsteps < mcp.thermalization_; ++thermsteps){
        std::vector<double> energies(betas.size());
        for (auto ii = 0u; ii < betas.size(); ++ii){
            mcs[ii].update_(sweeps);
            energies[ii] = mcs[ii].energy_();
        }
        pt.update_(energies);
        betas = pt.reorderedBetas_();
        
        for (auto ii = 0u; ii < betas.size(); ++ii)
            mcs[ii].beta_ = betas[ii];
    }
    #if 1
    std::cout << "Thermalized\n";
    for (auto meassteps = 0u; meassteps < mcp.measurement_; ++meassteps){
        std::vector<double> energies(betas.size());
        for (auto ii = 0u; ii < betas.size(); ++ii){
            mcs[ii].update_(sweeps);
            energies[ii] = mcs[ii].energy_();
        }
        os << pt.variableOrdered_(energies);
        pt.update_(energies);
        betas = pt.reorderedBetas_();
        for (auto ii = 0u; ii < betas.size(); ++ii)
            mcs[ii].beta_ = betas[ii];
    }
    #endif
}