#ifndef CLASSMAG_BASE_NUMERICS_HPP
#define CLASSMAG_BASE_NUMERICS_HPP

#include <cmath>
#include <vector>

#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    double pi();
    
    std::vector<double> shanks(std::vector<double> &x);
    
    template<unsigned int dimension>
    std::vector<geometry::Euclidean<dimension>> integerSweepPositive(unsigned int nmax){
        auto n_elems = static_cast<unsigned int>(pow(nmax,dimension + 1)) - 1u;
        auto result = std::vector<geometry::Euclidean<dimension>>(n_elems);
        for (auto ii = 0u; ii < n_elems; ++ii){
            for (auto jj = 0u; jj < dimension; ++jj){
                auto numerator = ii + 1u;
                auto denominator = static_cast<unsigned int>(pow(nmax,jj));
                result[ii][jj] = (numerator/denominator) % nmax;
            }
        }
        return result;
    }
    

    template<unsigned int dimension>
    std::vector<geometry::Euclidean<dimension>> integerSweepFull(unsigned int nmax){
        auto offsets = integerSweepPositive<dimension>(2u * nmax);
        auto corner = geometry::Euclidean<dimension>();
        corner.fill(-1.0 * nmax);
        auto result = std::vector<geometry::Euclidean<dimension>>(offsets.size() + 1);
        result[0] = corner;
        for (auto ii = 0u; ii < result.size(); ++ii){
            result[ii + 1] = corner + offsets[ii];
        }
        return result;
    }
    
}

#endif