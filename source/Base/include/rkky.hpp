#ifndef CLASSMAG_BASE_RKKY_HPP
#define CLASSMAG_BASE_RKKY_HPP

#include <math.h>
#include <functional>
#include <utility>

#include "Geometry/include/lattice.hpp"
#include "Base/include/spinStructure.hpp"

namespace classmag::base{
    std::function<double(const unsigned int, const unsigned int)>
        rkkyInteraction(
            const double kf,
            const geometry::Lattice<3> lattice,
            const double cutoff);
}

#endif