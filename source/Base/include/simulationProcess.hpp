#include "spinStructure.hpp"
#include "couplingLookup.hpp"

namespace classmag::base{
    
    template <unsigned int spinDimension>
    class SimulationProcess{
        public:
        SimulationProcess();
        void update_();
        void update_(const unsigned int n_times);
    };

}