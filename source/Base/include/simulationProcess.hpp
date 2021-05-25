#ifndef CLASSMAG_BASE_SIMULATIONPROCESS_HPP
#define CLASSMAG_BASE_SIMULATIONPROCESS_HPP

#include <functional>
#include <vector>
#include <stdexcept>

#include "spinStructure.hpp"
#include "couplingLookup.hpp"

namespace classmag::base{
    
    template<unsigned int spinDimension>
    class SimulationBase{
        public:
        SimulationBase(const SimulationBase &sb) = delete;

        virtual void update_(){

        }

        virtual std::vector<double> measure_() const{
            std::vector<double> result = {0.0};
            return result;
        }

        void update_(const unsigned int n_times){
            for (unsigned int ii = 0; ii < n_times; ++ii)
                update_();
        }

        void setSpinstructure_(SpinStructure<spinDimension>& spin){
            if (spin.size() != n_sites_){
                std::cerr << "Element number mismatch in spin initialization\n";
                abort();
            }

            spin_ = spin;

        }
        
        protected:
        SimulationBase(unsigned int n_sites):
        n_sites_(n_sites)
        {
            spin_.resize(n_sites);
        }

        SpinStructure<spinDimension> spin_;
        const unsigned int n_sites_;
    };
}

#endif