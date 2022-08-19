#ifndef CLASSMAG_MONTECARLO_MONTECARLOBASE_HPP
#define CLASSMAG_MONTECARLO_MONTECARLOBASE_HPP

#include <random>

#include "Base/include/simulationProcess.hpp"
#include "Base/include/orderParameter.hpp"

namespace classmag::montecarlo{
    template <unsigned int spinDimension>
    class MonteCarloBase : public base::SimulationBase<spinDimension>{
        public:
        void seed_(const int seed){
            rng_ = std::mt19937(seed);
        }

        void printSpin_(unsigned int site) const{
            for (unsigned int ii = 0; ii < spinDimension; ++ii)
                std::cout << this->spin_[site][ii] << " ";
            std::cout << "\n";
        }

        void printM_() const{
            for (unsigned int ii = 0; ii < spinDimension; ++ii){
                auto total = 0.0;
                for (auto s : this->spin_)
                    total += s[ii];
                std::cout << total << " ";
            }
            std::cout << "\n";
        }

        template<unsigned int latDim>
        void printSpins_(std::stringstream& stream, geometry::Lattice<latDim>&lattice) const{
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                for (auto s : this->spin_[ii])
                    stream << s << " ";
                for (auto x : lattice.position_(ii))
                    stream << x << " ";
                stream << "\n";
            }
        }

        void addOrderParameter_(
            const base::OrderParameter<spinDimension> &op){
            orderParameters_.push_back(op);
        }

        void addOrderParameter_(
            const base::OrderParameters<spinDimension> &op){
            for (auto o : op)
                addOrderParameter_(o);
        }

        void addOrderParameter_(
            const std::vector<base::OrderParameter<spinDimension>> &op){
            for (auto o : op)
                addOrderParameter_(o);
        }

        virtual std::vector<double> measure_() const {
            std::vector<double> results;
            for (auto orderParameter : orderParameters_)
                results.push_back(orderParameter.measure_(this->spin_));
            return results;
        }

        std::vector<std::string> measurementNames_() const
        {
            std::vector<std::string> results;
            for (auto measure : orderParameters_)
                results.push_back(measure.name_);
            return results;
        }

        double beta_ = 1.0;
        
        protected:
        MonteCarloBase(const unsigned int n):
        base::SimulationBase<spinDimension>(n){};

        std::mt19937 rng_ = std::mt19937(0);
        std::uniform_real_distribution<double> distr_ = 
            std::uniform_real_distribution<double>(0,1);

        private:
        base::OrderParameters<spinDimension> orderParameters_;
    };
}

#endif