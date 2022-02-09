#ifndef CLASSMAG_BASE_DEFAULT_ORDERPARAMETERS_HPP
#define CLASSMAG_BASE_DEFAULT_ORDERPARAMETERS_HPP

#include "orderParameter.hpp"

#include <vector>

namespace classmag::base
{
    OrderParameter<3> MagnetizationElement(unsigned int cartesian_index);

    std::vector<OrderParameter<3>> Magnetization();

    OrderParameter<3> BipartiteStaggeredMagnetization
    (
        const std::vector<unsigned int>& sublattice_partition_indices
    );
}

#endif