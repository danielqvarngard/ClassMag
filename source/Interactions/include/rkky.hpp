#ifndef CLASSMAG_MODELS_RKKY_HPP
#define CLASSMAG_MODELS_RKKY_HPP

#include <utility>
#include <math.h>

#include "Geometry/include/lattice.hpp"
#include "Interactions/include/spinStructure.hpp"

namespace classmag::models{
    std::function<double(const unsigned int, const unsigned int)>
        rkkyInteraction(
            const double kf,
            const geometry::Lattice<3> lattice,
            const double cutoff);
}

#endif