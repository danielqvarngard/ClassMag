#ifndef CLASSMAG_MODELS_SPINSTRUCTURE_HPP
#define CLASSMAG_MODELS_SPINSTRUCTURE_HPP

#include <vector>

namespace classmag::models{

    template<unsigned int spinDimension>
    class SpinStructure : 
        public std::vector<geometry::Euclidean<spinDimension>>{
        
    };
}


#endif