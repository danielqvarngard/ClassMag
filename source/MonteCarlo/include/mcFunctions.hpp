#ifndef CLASSMAG_MONTECARLO_MCFUNCTIONS_HPP
#define CLASSMAG_MONTECARLO_MCFUNCTIONS_HPP

#include <cmath>

#include "Geometry/include/euclidean.hpp"


namespace classmag::montecarlo{
    double exponential(double x);
    double logarithm(double x);
    double boltzmannFactor(double beta, double delta_E);
    double heatBathCosine(double fieldStrength, double z);
}

#endif