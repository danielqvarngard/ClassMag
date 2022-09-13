#ifndef CLASSMAG_BASE_DIPOLE_SPHERICAL
#define CLASSMAG_BASE_DIPOLE_SPHERICAL

#include "dipole_real_space.hpp"
#include "dipole_reciprocal_space.hpp"
#include "dipole_profile.hpp"

namespace classmag::base
{

    geometry::Matrix<3,3> dipole_spherical_matrix
    (
        const unsigned int site1,
        const unsigned int site2,
        const DipoleProfile& ep
    );

     geometry::Matrix<3,3> dipole_fishy_matrix
    (
        const unsigned int site1,
        const unsigned int site2,
        const DipoleProfile& ep
    );
}

#endif