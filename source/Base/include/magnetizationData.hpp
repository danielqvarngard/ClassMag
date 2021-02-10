#ifndef CLASSMAG_BASE_MAGNETIZATIONDATA_HPP
#define CLASSMAG_BASE_MAGNETIZATIONDATA_HPP

#include <vector>
#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    template <unsigned int spinDimension>
    class MagnetizationData : public std::vector<geometry::Euclidean<spinDimension>>
    {
        public:
        void fill_(const double x){
            for (auto ii = 0u; ii < (*this).size(); ++ii)
                (*this)[ii].fill_(x);
        }
    };
}


#endif