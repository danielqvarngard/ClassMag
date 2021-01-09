#include <vector>
#include <random>

#include "mcFunctions.hpp"

namespace classmag::montecarlo{
    struct PermutationEntry{
        unsigned int variableIndex;
        unsigned int permutationIndex;
    };

    class PermutationManager{
        private:
        std::vector<PermutationEntry> sn_;

        public:
        PermutationManager(unsigned int size);

        void switch_(unsigned int a, unsigned int b);

        unsigned int variable_(unsigned int processIndex);
        unsigned int process_(unsigned int variableIndex);
    };

    class ParallelTemperer{
        private:
        PermutationManager pm_;
        std::vector<double> temperatures_;

        std::mt19937 mt_ = std::mt19937(0);
        std::uniform_real_distribution<double> distr_ =
            std::uniform_real_distribution<double>(0.0,1.0);
        

        public:
        ParallelTemperer(const std::vector<double> &temperatures);
        
        virtual void update_(const std::vector<double> &energies);
    };
}