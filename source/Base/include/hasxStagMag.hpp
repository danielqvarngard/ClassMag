#ifndef CLASSMAG_BASE_HASXSTAGMAG_HPP
#define CLASSMAG_BASE_HASXSTAGMAG_HPP

#include "orderParameter.hpp"

namespace classmag::base{
    std::pair<OrderParameter<3>, OrderParameter<3>> chasxStagMag_split(
        const geometry::Lattice<3> &lattice);

    std::pair<OrderParameter<3>, OrderParameter<3>> chasxMag_split(
        const geometry::Lattice<3> &lattice);
}

#endif