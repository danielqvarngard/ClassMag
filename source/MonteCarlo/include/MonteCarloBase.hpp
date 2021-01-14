#ifndef CLASSMAG_MONTECARLO_MONTECARLOBASE_HPP
#define CLASSMAG_MONTECARLO_MONTECARLOBASE_HPP

#include <random>

#include "Base/include/simulationProcess.hpp"

namespace classmag::montecarlo{
    template <unsigned int spinDimension>
    class MonteCarloBase : public base::SimulationBase<spinDimension>{
        public:
        double beta_ = 1.0;
        private:
        std::mt19937 rng_ = std::mt19937(0);
        std::uniform_real_distribution<double> distr_ = 
            std::uniform_real_distribution<double>(0,1);
    };
}

#endif