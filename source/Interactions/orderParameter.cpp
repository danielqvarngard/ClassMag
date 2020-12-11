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
}