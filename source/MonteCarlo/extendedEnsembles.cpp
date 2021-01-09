#include "include/extendedEnsembles.hpp"

namespace classmag::montecarlo{
    PermutationManager::PermutationManager(unsigned int size):
    sn_(std::vector<PermutationEntry>(size))
    {   
        auto p = 0u;
        for(auto e : sn_){
            e.permutationIndex = p;
            e.variableIndex = p;
            ++p;
        }
    }

    void PermutationManager::switch_(unsigned int a, unsigned int b){
        auto temp = sn_[a].variableIndex;

        sn_[sn_[a].variableIndex].permutationIndex = b;
        sn_[sn_[b].variableIndex].permutationIndex = a;
        sn_[a].variableIndex = sn_[b].variableIndex;
        sn_[b].variableIndex = temp;
    }

    unsigned int PermutationManager::process_(unsigned int variableIndex){
        return sn_[variableIndex].permutationIndex;
    }

    unsigned int PermutationManager::variable_(unsigned int processIndex){
        return sn_[processIndex].variableIndex;
    }

    ParallelTemperer::ParallelTemperer(const std::vector<double> &temperatures):
    temperatures_(temperatures),
    pm_(PermutationManager(temperatures.size()))
    {

    }

    void ParallelTemperer::update_(const std::vector<double> &energies){
        for (auto ii = 0u; ii < temperatures_.size() - 1; ii = ii + 2){
            auto deltaBeta = 1.0/temperatures_[ii+1] - 1.0/temperatures_[ii];
            auto deltaE = energies[pm_.process_(ii+1)] - energies[pm_.process_(ii)];
            if (boltzmannFactor(deltaBeta,deltaE) > distr_(mt_))
                pm_.switch_(ii + 1, ii);
        }

        for (auto ii = 1u; ii < temperatures_.size() - 1; ii = ii + 2){
            pm_.switch_(ii, ii + 1);
        }
    }

}