#include "include/orderParameter.hpp"
#include "Geometry/include/predefLattices.hpp"

namespace classmag::models{
    std::pair<geometry::Lattice<3>, OrderParameter<3>> DefaultTsaiSim(
        const std::array<unsigned int,3> &systemSize){

        auto lattice = geometry::has0(systemSize);
        auto orderParameter = clusterOrder<3,3>(lattice);

        std::pair<geometry::Lattice<3>, OrderParameter<3>> result =
            {lattice, orderParameter};
        return result;
    };

    std::pair<OrderParameter<3>, OrderParameter<3>> has100Params(
        const std::array<unsigned int,3> &systemSize){
            auto lattice = geometry::has100(systemSize);
            std::function<double(const SpinStructure<3>&)> 
                shellOrder = 
                [lattice](const SpinStructure<3> &spins){
                auto result = 0.0;
                for (unsigned int site = 0; site < spins.size(); ++site){
                    if ((site + 1) % 13 != 0){
                        auto loopIndices = lattice.neighborCellSites_(site);
                        auto n_neighbors = static_cast<double>(loopIndices.size());
                        for (auto neighborSite : loopIndices)
                            result += spins[site] * spins[neighborSite]/n_neighbors;
                    }
                }
                return result;
            };

            std::function<double(const SpinStructure<3>&)> 
                centerOrder = 
                [lattice](const SpinStructure<3> &spins){
                auto result = 0.0;
                for (unsigned int site = 12; site < spins.size(); site = site + 13){
                    auto loopIndices = lattice.neighborCellSites_(site);
                    auto n_neighbors = static_cast<double>(loopIndices.size());
                    for (auto neighborSite : loopIndices)
                        result += spins[site] * spins[neighborSite]/n_neighbors;
                }
                return result;
            };

            auto opShell = OrderParameter<3>("Shell order", shellOrder);
            auto opCenter = OrderParameter<3>("Center order", centerOrder);
            auto opPair = std::pair<OrderParameter<3>,OrderParameter<3>>({opShell,opCenter});
            return opPair;
        };
}