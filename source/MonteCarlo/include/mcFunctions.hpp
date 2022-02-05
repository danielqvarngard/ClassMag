#ifndef CLASSMAG_MONTECARLO_MCFUNCTIONS_HPP
#define CLASSMAG_MONTECARLO_MCFUNCTIONS_HPP

#include <cmath>
#include <vector>

#include "Geometry/include/euclidean.hpp"


namespace classmag::montecarlo{
    double squareRoot(double x);
    double pi();
    double exponential(double x);
    double logarithm(double x);
    double boltzmannFactor(double beta, double delta_E);
    double heatBathCosine(double fieldStrength, double z);
    std::vector<double> invert_vector_elements(const std::vector<double>& v);
}

#endif