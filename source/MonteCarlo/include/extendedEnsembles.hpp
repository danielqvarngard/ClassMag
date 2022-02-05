#ifndef CLASSMAG_MONTECARLO_EXTENDEDENSEMBLES_HPP
#define CLASSMAG_MONTECARLO_EXTENDEDENSEMBLES_HPP

#include <vector>
#include <random>

#include "mcFunctions.hpp"
#include "Parallelism/include/ArrayMessage.hpp"

namespace classmag::montecarlo{
    struct PermutationEntry{
        unsigned int variableIndex;
        unsigned int processIndex;
    };

    class PermutationManager{
        public:
        PermutationManager(const PermutationManager &) = delete;

        inline unsigned int variable_(const unsigned int processIndex) const;
        std::vector<unsigned int> variable_();
        inline unsigned int process_(const unsigned int variableIndex) const;

        std::vector<double> processOrdered_(const std::vector<double> &variableordered);
        std::vector<double> variableOrdered_(const std::vector<double> &processordered);
        parallelism::ArrayMessage variableOrdered_(const parallelism::ArrayMessage &processordered);
        protected:
        PermutationManager(unsigned int size);
        inline void switchProcess_(unsigned int a, unsigned int b);

        private:
        std::vector<PermutationEntry> sn_;

    };

    class ParallelTemperer : public PermutationManager{
        private:
        void seed_(const int seed);
        std::vector<double> betas_;
        std::vector<double> acceptance_rates_;

        std::mt19937 mt_ = std::mt19937(0);
        std::uniform_real_distribution<double> distr_ =
            std::uniform_real_distribution<double>(0.0,1.0);
        
        public:
        ParallelTemperer(const std::vector<double> &betas);
        
        virtual void update_(const std::vector<double> &energies);

        std::vector<double> reorderedBetas_();
        std::vector<double> get_acceptance_rates_();
    };
}

#endif