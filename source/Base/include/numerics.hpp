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
        auto result = integerSweepPositive<dimension>(nmax);
        auto negative = result;
        for (auto ii = 0u; ii < negative.size(); ++ii){
            negative[ii] *= -1.0;
        }
        result.insert(std::end(result), std::begin(negative), std::end(negative));
        auto origin = geometry::Euclidean<dimension>();
        origin.fill(0.0);
        result.push_back(origin);
        return result;
    }
    
}

#endif