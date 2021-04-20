#include <array>
#include <iostream>
#include <algorithm>
#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/mcProfile.hpp"

using namespace classmag;

template<unsigned int dimension>
void printLattice(const geometry::Lattice<dimension> lattice){
    for (unsigned int ii = 0; ii < lattice.n_sites_(); ++ii){
        auto e = lattice.position_(ii);
        for (auto d : e)
            std::cout << d << " ";
        std::cout << "\n";
    }
};

int main(){

    

    #if 0
    montecarlo::VectorModel_Profile mcp;
    auto n_thermalize = 10000;
    auto n_overrelax = 1;
    auto n_measure = 10000;
    auto n_skip = 1;
    auto n_resamples = 100;
    auto L = 2u;
    const auto systemSize = std::array<unsigned int, 3>({L,L,L});

    /* bcc lattice for testing: */
    
    auto sublattice1 = geometry::cubicLattice<3>(systemSize);
    auto lattice = geometry::Lattice<3>(sublattice1);
    auto sublattice2 = geometry::cubicLattice<3>(systemSize);
    auto e = geometry::Euclidean<3>({0.5, 0.5, 0.5});
    sublattice2.decorate_({e});

    lattice.append_(sublattice2);
    const auto interaction = base::nearestNeighbor(-1.0,lattice,1.01);
    
    
    
    auto lattice = geometry::chas100(systemSize);
    const auto interaction = base::nearestNeighbor(-1.0,lattice,0.385);
   
    mcp.measurement_ = n_measure;
    mcp.thermalization_ = n_thermalize;
    mcp.skips_ = n_skip;
    mcp.n_sites_ = lattice.n_sites_();
    mcp.overrelax_ = n_overrelax;
    mcp.seed_ = 137;
    mcp.getPartitions_(lattice);

    std::cout << mcp.n_sites_ << "\n";

    std::string dir = "../out/";
    std::string filename = dir + "testMC_nanbug_";
    auto seed = 0;
    auto mc = montecarlo::VectorModelManager<3>(
        mcp,
        interaction);

    mc.printCoordinations_();

    printLattice(lattice);

    std::cout << lattice.n_decorations_() << "\n";
    #endif
};