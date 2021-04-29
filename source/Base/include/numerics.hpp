#ifndef CLASSMAG_BASE_NUMERICS_HPP
#define CLASSMAG_BASE_NUMERICS_HPP

#include <cmath>
#include <vector>

#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    double pi();
    
    std::vector<double> shanks(std::vector<double> &x);

    std::vector<double> linspace(
        const double min, 
        const double max, 
        const unsigned int stepCount);
    
    template<unsigned int dimension>
    std::vector<geometry::Euclidean<dimension>> integerSweepPositive(unsigned int nmax){
        auto n_elems = static_cast<unsigned int>(pow(nmax + 1, dimension)) - 1u;
        auto result = std::vector<geometry::Euclidean<dimension>>(n_elems);
        for (auto ii = 0u; ii < n_elems; ++ii){
            for (auto jj = 0u; jj < dimension; ++jj){
                auto numerator = ii + 1u;
                auto denominator = static_cast<unsigned int>(pow(nmax + 1,jj));
                result[ii][jj] = (numerator/denominator) % (nmax + 1);
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

    template<unsigned int dimension>
    std::vector<geometry::Euclidean<dimension>> integerSweepExclude(unsigned int nmax){
        auto offsets = integerSweepPositive<dimension>(2u * nmax);
        auto corner = geometry::Euclidean<dimension>();
        corner.fill(-1.0 * nmax);
        auto result = std::vector<geometry::Euclidean<dimension>>(offsets.size());
        result[0] = corner;
        auto originIndex = 0u;
        for (auto ii = 0u; ii < dimension; ++ii){
            originIndex += static_cast<unsigned int>(pow(2u*nmax + 1u,ii)) * nmax;
        }
        originIndex -= 1u;
        for (auto ii = 0u; ii < originIndex; ++ii){
            result[ii + 1] = corner + offsets[ii];
        }
        for (auto ii = originIndex + 1u; ii < offsets.size(); ++ii){
            result[ii] = corner + offsets[ii];
        }
        return result;
    }
    
}

#endif