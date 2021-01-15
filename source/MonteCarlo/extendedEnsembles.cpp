#include "include/extendedEnsembles.hpp"

namespace classmag::montecarlo{
    PermutationManager::PermutationManager(unsigned int size):
    sn_(std::vector<PermutationEntry>(size))
    {   
        auto p = 0u;
        for(auto e : sn_){
            e.processIndex = p;
            e.variableIndex = p;
            ++p;
        }
    }

    inline void PermutationManager::switchProcess_(unsigned int var_a, unsigned int var_b){
        auto temp = sn_[var_a].processIndex;

        sn_[sn_[var_a].processIndex].variableIndex = var_b;
        sn_[sn_[var_b].processIndex].variableIndex = var_a;
        sn_[var_a].processIndex = sn_[var_b].processIndex;
        sn_[var_b].processIndex = temp;
    }

    unsigned int PermutationManager::process_(unsigned int variableIndex){
        return sn_[variableIndex].processIndex;
    }

    unsigned int PermutationManager::variable_(unsigned int processIndex){
        return sn_[processIndex].variableIndex;
    }

    ParallelTemperer::ParallelTemperer(const std::vector<double> &betas):
    PermutationManager{betas.size()},
    betas_(betas)
    {

    }

    void ParallelTemperer::update_(const std::vector<double> &energies){
        //TODO: add size check / refactor as template 
        for (auto ii = 0u; ii < betas_.size() - 1; ii = ii + 2){
            auto deltaBeta = betas_[ii+1] - betas_[ii];
            auto deltaE = energies[process_(ii+1)] - energies[process_(ii)];
            if (boltzmannFactor(deltaBeta,deltaE) > distr_(mt_))
                switchProcess_(ii + 1, ii);
        }

        for (auto ii = 1u; ii < betas_.size() - 1; ii = ii + 2){
            auto deltaBeta = betas_[ii+1] - betas_[ii];
            auto deltaE = energies[process_(ii+1)] - energies[process_(ii)];
            if (boltzmannFactor(deltaBeta,deltaE) > distr_(mt_))
                switchProcess_(ii + 1, ii);
        }
    }

}