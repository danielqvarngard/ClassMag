#include "ptProcs.hpp"

namespace classmag::environments{
    int mpiPT_hub(
        const std::vector<double> &betas, 
        const unsigned int orderCount,
        const unsigned int measurementCount){
        if (orderCount == 0)
            return mpiPT_hub(betas);

        parallelism::Hub msg;
        parallelism::VectorTarget vt(orderCount);
        montecarlo::ParallelTemperer ptm(betas);
        
        msg.scatterDoubles_(betas,betaChannel);

        for (auto ii = 0u; ii < measurementCount; ++ii){
            std::vector<double> energies;
            msg.gatherDoubles_(energies, energyChannel);
            
        }



        return 0;
    };
}