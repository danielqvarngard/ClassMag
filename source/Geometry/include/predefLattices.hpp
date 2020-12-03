#ifndef CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP
#define CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP

#include "lattice.hpp"

namespace classmag::geometry{
    std::vector<Euclidean<3>> icosahedralCluster();

    Lattice<2> squareLattice(const std::array<unsigned int,2> &systemSize);
}

#endif