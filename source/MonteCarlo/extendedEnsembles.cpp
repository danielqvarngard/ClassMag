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
        
    }

}