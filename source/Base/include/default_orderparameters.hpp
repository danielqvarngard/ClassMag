#ifndef CLASSMAG_BASE_DEFAULT_ORDERPARAMETERS_HPP
#define CLASSMAG_BASE_DEFAULT_ORDERPARAMETERS_HPP

#include "orderParameter.hpp"

#include <vector>

namespace classmag::base
{
    OrderParameter<3> MagnetizationElements(unsigned int cartesian_index);
    std::vector<OrderParameter<3>> Magnetization
    (

    );

    OrderParameter<3> BipartiteStaggeredMagnetization
    (
        unsigned int cartesian_index,
        std::vector<unsigned int> sublattice_partition_indices
    );
    
    std::vector<OrderParameter<3>> BipartiteStaggeredMagnetization
    (
        std::vector<unsigned int> sublattice_partition_indices
    );
}

#endif