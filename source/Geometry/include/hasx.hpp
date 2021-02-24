#ifndef CLASSMAG_GEOMETRY_HASX_HPP
#define CLASSMAG_GEOMETRY_HASX_HPP

#include "disorderedLattices.hpp"

namespace classmag::geometry{
    Lattice<3> hasx(const std::array<unsigned int, 3> &size, const double occupancy){
        auto result = chas0_bipartite(size);
        auto disp = Euclidean<3>({0.5, 0.5, 0.5});

        result.append_(disorderedCubic(size, occupancy));
        result.append_(disorderedCubic(size, occupancy, disp));
        return result;
    }
}

#endif