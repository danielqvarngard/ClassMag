#ifndef CLASSMAG_GEOMETRY_RECIPROCAL_HPP
#define CLASSMAG_GEOMETRY_RECIPROCAL_HPP

#include "euclidean.hpp"

namespace classmag::geometry{
    std::array<Euclidean<3>,3> reciprocalBasis(const std::array<Euclidean<3>,3> &bravais);
}

#endif