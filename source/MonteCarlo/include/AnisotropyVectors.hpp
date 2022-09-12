#ifndef CLASSMAG_ANISOTROPYVECTORS_HPP
#define CLASSMAG_ANISOTROPYVECTORS_HPP

#include <vector>
#include <array>

#include "Geometry/include/euclidean.hpp"

namespace classmag::montecarlo{
    template<unsigned int spinDim>
    class EasyPlaneVectors : 
        public std::vector<
            std::array<geometry::Euclidean<spinDim>,2>>{
    };

    template<unsigned int spinDim>
    class EasyAxisVectors :
        public std::vector<
            std::vector<
                geometry::Euclidean<spinDim>>>{
    };

}

#endif