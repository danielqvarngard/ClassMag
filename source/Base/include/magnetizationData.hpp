#ifndef CLASSMAG_BASE_MAGNETIZATIONDATA_HPP
#define CLASSMAG_BASE_MAGNETIZATIONDATA_HPP

#include <array>
#include <iostream>
#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    template <unsigned int spinDimension, unsigned int subLattices>
    class MagnetizationData : 
    public std::array<geometry::Euclidean<spinDimension>, subLattices>
    {
        public:
        void fill_(const double x){
            for (auto ii = 0u; ii < (*this).size(); ++ii)
                (*this)[ii].fill(x);
        }

        MagnetizationData<spinDimension, subLattices>& operator=(
            const MagnetizationData<spinDimension,subLattices> &x){
            for (auto ii = 0u; ii < subLattices; ++ii){
                (*this)[ii] = x[ii];
            }
            return *this;
        }

        MagnetizationData<spinDimension, subLattices>& operator=(
            const std::vector<geometry::Euclidean<spinDimension>> &x){
            for (auto ii = 0u; ii < subLattices; ++ii){
                (*this)[ii] = x[ii];
            }
            return *this;
        }

        MagnetizationData<spinDimension, subLattices>& operator+=(
            const MagnetizationData<spinDimension,subLattices> &x){
            for (auto ii = 0u; ii < subLattices; ++ii){
                (*this)[ii] += x[ii];
            }
            return *this;
        }

        MagnetizationData<spinDimension, subLattices>& operator*=(
            const double &x){
            for (auto ii = 0u; ii < subLattices; ++ii){
                (*this)[ii] *= x;
            }
            return *this;
        }
    };

    template<unsigned int spinDimension, unsigned int subLattices> 
    MagnetizationData<spinDimension, subLattices> operator*(
        const double d, 
        const MagnetizationData<spinDimension, subLattices>& src){
        MagnetizationData<spinDimension, subLattices> result;
        result = src;
        result *= d;
        return result;
    }


    template<unsigned int spinDimension, unsigned int subLattices>
    MagnetizationData<spinDimension, subLattices> operator-(
        const MagnetizationData<spinDimension, subLattices>& src1,
        const MagnetizationData<spinDimension, subLattices>& src2){
        
        MagnetizationData<spinDimension, subLattices> result;
        result = src1; 
        result += (-1.0)*src2;
        return result;
    }

    template<unsigned int spinDimension, unsigned int subLattices>
    MagnetizationData<spinDimension, subLattices> operator*(
        const MagnetizationData<spinDimension, subLattices>& src1,
        const MagnetizationData<spinDimension, subLattices>& src2){
        
        MagnetizationData<spinDimension, subLattices> result;
        result = src1; 
        for (auto ii = 0u; ii < subLattices; ++ii){
            for (auto jj = 0u; jj < spinDimension; ++jj)
                result[ii][jj] *= src2[ii][jj];
        }
        return result;
    }
}


#endif