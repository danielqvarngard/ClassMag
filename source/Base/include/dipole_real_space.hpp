#ifndef CLASSMAG_BASE_DIPOLE_REAL_SPACE_HPP
#define CLASSMAG_BASE_DIPOLE_REAL_SPACE_HPP

#include "Geometry/include/matrix.hpp"
#include "numerics.hpp"

namespace classmag::base
{
    geometry::Matrix<3,3> dipole_real_space_sum
    (
        const geometry::Euclidean<3> displacement,
        const double alpha,
        const unsigned int mirrors
    );

    geometry::Matrix<3,3> dipole_real_space_term
    (
        const geometry::Euclidean<3> folded_displacement,
        const double alpha
    );

    double dipole_real_scalar
    (
        const double folded_distance, 
        const double alpha
    );

    double dipole_real_matrix_common_factor
    (
        const double folded_distance, 
        const double alpha
    );

    geometry::Matrix<3,3> dipole_real_ext_prod
    (
        const geometry::Euclidean<3> folded_displacement
    );
}

#endif