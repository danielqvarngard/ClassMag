#ifndef CLASSMAG_BASE_NUMERICS_HPP
#define CLASSMAG_BASE_NUMERICS_HPP

#include <cmath>
#include <vector>

#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    inline double pi()
    {
        return M_PI;
    };
    
    std::vector<double> shanks(std::vector<double> &x);

    std::vector<double> linspace(
        const double min, 
        const double max, 
        const unsigned int stepCount
    );

    std::vector<double> linspace_open(
        const double min, 
        const double max, 
        const unsigned int stepCount
    );

    std::vector<double> logspace(
        const double min,
        const double max,
        const unsigned int step_count
    );

    double trapz_term(const double d0, const double d1, const double dx);
    double euler_forward(const double d0, const double d1, const double dx);
    
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
        for (auto ii = 0u; ii < offsets.size(); ++ii){
            result[ii + 1] = corner + offsets[ii];
        }
        return result;
    }

    template<unsigned int dimension>
    std::vector<geometry::Euclidean<dimension>> integerSweepExclude(const unsigned int nmax){
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
    
    inline unsigned int maximum_z3plus_index(const unsigned int maximum_entry)
    {
        unsigned int result = 
            maximum_entry * (maximum_entry + 1) * (maximum_entry + 1) + 
            maximum_entry * (maximum_entry + 1) +
            maximum_entry + 1;

        return result;
    }

    inline std::array<unsigned int, 3> z3plus_entry(const unsigned int index, const unsigned int max)
    {
        auto n1 = index % (max + 1);
        auto n2 = (index/(max + 1)) % (max + 1);
        auto n3 = (index/((max + 1u) * (max + 1u))) % (max + 1);
        
        auto result = std::array<unsigned int, 3>({n1, n2, n3});
        return result;
    }
}

#endif