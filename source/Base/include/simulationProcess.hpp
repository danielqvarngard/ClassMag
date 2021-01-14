#ifndef CLASSMAG_BASE_SIMULATIONPROCESS_HPP
#define CLASSMAG_BASE_SIMULATIONPROCESS_HPP

#include <functional>

#include "spinStructure.hpp"
#include "couplingLookup.hpp"

namespace classmag::base{
    
    template <unsigned int spinDimension>
    class SimulationBase{
        public:
        SimulationBase(){
        };

        virtual void update_(){

        };

        void update_(const unsigned int n_times){
            for (unsigned int ii = 0; ii < n_times; ++ii)
                update_();
        };
    };

}

#endif