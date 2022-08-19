#include "include/dipole_spherical.hpp"

namespace classmag::base
{
    geometry::Matrix<3,3> dipole_spherical_matrix
    (
        const unsigned int site1,
        const unsigned int site2,
        const DipoleProfile& ep
    ){
        auto r = ep.lattice_.position_(site1) - ep.lattice_.position_(site2);
        r = 1.0/ep.length_ * r;
        auto real_space_component = dipole_real_space_sum(r, ep.alpha_, ep.realMirrors_);
        auto rec_space_component = 
            dipole_reciprocal_space_sum(r, ep.alpha_, ep.recMirrors_);

        auto result = real_space_component + rec_space_component;
        return (-1.0 * ep.magnitude_/pow(ep.length_, 3.0)) * result;
    };
}