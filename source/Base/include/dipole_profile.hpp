#ifndef CLASSMAG_BASE_DIPOLE_PROFILE_HPP
#define CLASSMAG_BASE_DIPOLE_PROFILE_HPP

#include "Geometry/include/lattice.hpp"

namespace classmag::base 
{
    struct DipoleProfile{
        public:
        DipoleProfile(geometry::Lattice<3> lat):
            lattice_(lat)
            {

        };
        geometry::Lattice<3> lattice_;
        double alpha_ = 1.0;
        double magnitude_ = 1.0;
        double length_;
        unsigned int realMirrors_ = 0u;
        unsigned int recMirrors_ = 0u;
    };
}

#endif