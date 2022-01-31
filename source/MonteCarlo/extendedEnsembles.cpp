#include "include/extendedEnsembles.hpp"

namespace classmag::montecarlo{
    PermutationManager::PermutationManager(unsigned int size):
    sn_(std::vector<PermutationEntry>(size))
    {   
        for(auto p = 0u; p < sn_.size(); ++p){
            sn_[p].processIndex = p;
            sn_[p].variableIndex = p;
        }
    }

    inline void PermutationManager::switchProcess_(unsigned int var_a, unsigned int var_b){
        auto temp = sn_[var_a].processIndex;

        sn_[sn_[var_a].processIndex].variableIndex = var_b;
        sn_[sn_[var_b].processIndex].variableIndex = var_a;
        sn_[var_a].processIndex = sn_[var_b].processIndex;
        sn_[var_b].processIndex = temp;
    }

    inline unsigned int PermutationManager::process_(const unsigned int variableIndex) const{
        return sn_[variableIndex].processIndex;
    }

    inline unsigned int PermutationManager::variable_(const unsigned int processIndex) const{
        return sn_[processIndex].variableIndex;
    }

    std::vector<unsigned int> PermutationManager::variable_(){
        std::vector<unsigned int> result(sn_.size());

        for (auto ii = 0u; ii < sn_.size(); ++ii){
            result[ii] = variable_(ii);
        }

        return result;
    }

    std::vector<double> PermutationManager::processOrdered_(
        const std::vector<double> &variableordered){
        std::vector<double> result(sn_.size());
        for (auto ii = 0u; ii < sn_.size(); ++ii){
            result[process_(ii)] = variableordered[ii];
        }
        return result;
    }

    std::vector<double> PermutationManager::variableOrdered_(
        const std::vector<double> &processordered){
        std::vector<double> result(sn_.size());
        for (auto ii = 0u; ii < sn_.size(); ++ii){
            result[variable_(ii)] = processordered[ii];
        }
        return result;
    }

    fileio::VectorTarget PermutationManager::variableOrdered_(
        const fileio::VectorTarget &processordered){
        fileio::VectorTarget result(processordered.messageLength_);
        for (auto ii = 0u; ii < processordered.messageLength_; ++ii){
            for (auto jj = 0u; jj < processordered.data_[ii].size(); ++jj)
                result.data_[variable_(ii)][jj] = processordered.data_[ii][jj];
        }
        return result;
    }

    ParallelTemperer::ParallelTemperer(const std::vector<double> &betas):
    PermutationManager(static_cast<unsigned int>(betas.size())),
    betas_(betas)
    {

    }

    void ParallelTemperer::seed_(const int seed){
        mt_ = std::mt19937(seed);
    }

    void ParallelTemperer::update_(const std::vector<double> &energies){
        //TODO: add size check / refactor as template 
        for (auto ii = 0u; ii < betas_.size() - 1; ii += 2){
            auto deltaBeta = betas_[ii+1] - betas_[ii];
            auto deltaE = energies[process_(ii+1)] - energies[process_(ii)];
            if (boltzmannFactor(-deltaBeta,deltaE) > distr_(mt_))
                switchProcess_(ii + 1, ii);
        }

        for (auto ii = 1u; ii < betas_.size() - 1; ii += 2){
            auto deltaBeta = betas_[ii+1] - betas_[ii];
            auto deltaE = energies[process_(ii+1)] - energies[process_(ii)];
            if (boltzmannFactor(-deltaBeta,deltaE) > distr_(mt_))
                switchProcess_(ii + 1, ii);
        }
    }

    std::vector<double> ParallelTemperer::reorderedBetas_(){
        auto permIndices = variable_();
        std::vector<double> result(betas_.size());
        for (auto ii = 0u; ii < betas_.size(); ++ii)
            result[ii] = betas_[permIndices[ii]];
        return result;
    }
}