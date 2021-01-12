#include <functional>

#include "spinStructure.hpp"
#include "couplingLookup.hpp"

namespace classmag::base{
    
    template <unsigned int spinDimension>
    class SimulationProcess{
        public:
        SimulationProcess(
            const unsigned int n_sites,
            const std::function<double(unsigned int, unsigned int)> interaction
            ):
            n_sites_(n_sites),
            lookup_(CouplingLookup(interaction)){
            
            spin_.resize(n_sites_);
        };

        virtual void update_(){

        };

        void update_(const unsigned int n_times){
            for (unsigned int ii = 0; ii < n_times; ++ii)
                update_();
        };

        private:
        const unsigned int n_sites_;
        SpinStructure<spinDimension> spin_;
        const CouplingLookup lookup_;
    };

}