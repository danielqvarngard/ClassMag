#include "ptProcs.hpp"

namespace classmag::environments{

    void MPI_PT_Hub_Thermalize
    (
        const PT_Input& input, 
        montecarlo::ParallelTemperer& ptm
    ){
        parallelism::CentralNodeMessenger msg; 
        for (auto ii = 0u; ii < input.thermalization_runs_; ++ii)
        {
            msg.scatterDoubles_(ptm.processOrdered_(input.betas_), betaChannel);
            std::vector<double> energies;
            msg.gatherDoubles_(energies,energyChannel);
            ptm.update_(energies);
        }    
    }

    PT_Return MPI_PT_Hub(const PT_Input& input){
        PT_Return result;
        montecarlo::ParallelTemperer ptm(input.betas_);
        MPI_PT_Hub_Thermalize(input, ptm);

        parallelism::CentralNodeMessenger msg;
        for (auto ii = 0u; ii < input.measurement_runs_; ++ii){
            std::vector<double> energies;
            msg.gatherDoubles_(energies, energyChannel); 
            ptm.update_(energies);
            msg.scatterDoubles_(ptm.processOrdered_(input.betas_), betaChannel);
        }
        return result;
    };

    PT_Return MPI_PT_Hub(
        const PT_Input& input,
        const unsigned int orderCount){
        if (orderCount == 0)
            return MPI_PT_Hub(input);
        
        montecarlo::ParallelTemperer ptm(input.betas_);
        parallelism::CentralNodeMessenger msg;
        parallelism::ArrayMessage vt(orderCount);
        PT_Return result;

        MPI_PT_Hub_Thermalize(input, ptm);
        ptm.reset_acceptance_rates();

        for (auto ii = 0u; ii < input.measurement_runs_; ++ii){
            msg.scatterDoubles_(ptm.processOrdered_(input.betas_), betaChannel);
            std::vector<double> energies;
            msg.gatherDoubles_(energies, energyChannel);
            msg.gatherDoubles_(vt,orderChannel);
            result.microstate_energies.push_back(ptm.variableOrdered_(energies));
            auto reordered = ptm.variableOrdered_(vt);
            result.microstate_variables.push_back(reordered);
            ptm.update_(energies);
        }

        auto accrate = ptm.get_acceptance_rates_();
        auto n = static_cast<double>(input.measurement_runs_);
        for (auto ii = 0; ii < accrate.size(); ++ii)
            accrate[ii] /= n;

        result.acceptance_rates = std::move(accrate);
        return std::move(result);
    };
}