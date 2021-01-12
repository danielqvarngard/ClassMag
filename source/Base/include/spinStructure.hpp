#ifndef CLASSMAG_BASE_SPINSTRUCTURE_HPP
#define CLASSMAG_BASE_SPINSTRUCTURE_HPP

#include <vector>

namespace classmag::base{

    template<unsigned int spinDimension>
    class SpinStructure : 
        public std::vector<geometry::Euclidean<spinDimension>>{
        
    };
}


#endif