#ifndef CLASSMAG_BASE_NUMERICS_HPP
#define CLASSMAG_BASE_NUMERICS_HPP

#include <cmath>
#include <vector>

#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    double pi();
    std::vector<double> shanks(std::vector<double> &x);
}

#endif