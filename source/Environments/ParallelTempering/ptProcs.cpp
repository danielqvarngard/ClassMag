#include "ptProcs.hpp"

namespace classmag::environments{
    int mpiPT_hub(
        const std::vector<double> &betas,
        const unsigned int measurementCount){
        
        auto fn = "out/" + fileio::datestamp();
        fileio::OStreamManager mcenFile(fn + ".mcen");
        fileio::OStreamManager mcopFile(fn + ".mcop");

        parallelism::Hub msg;
        montecarlo::ParallelTemperer ptm(betas);
        msg.scatterDoubles_(betas,betaChannel);

        for (auto ii = 0u; ii < measurementCount; ++ii){
            std::vector<double> energies;
            msg.gatherDoubles_(energies, energyChannel);
            ptm.update_(energies);
            msg.scatterDoubles_(ptm.processOrdered_(betas), betaChannel);
            mcenFile << ptm.variableOrdered_(energies);
        }

        return 0;
    };

    int mpiPT_hub(
        const std::vector<double> &betas, 
        const unsigned int measurementCount,
        const unsigned int orderCount){
        if (orderCount == 0)
            return mpiPT_hub(betas, measurementCount);

        auto fn = "out/" + fileio::datestamp();
        fileio::OStreamManager mcenFile(fn + ".mcen");
        fileio::OStreamManager mcopFile(fn + ".mcop");

        parallelism::Hub msg;
        fileio::VectorTarget vt(orderCount);
        montecarlo::ParallelTemperer ptm(betas);
        
        msg.scatterDoubles_(betas,betaChannel);

        for (auto ii = 0u; ii < measurementCount; ++ii){
            std::vector<double> energies;
            msg.gatherDoubles_(energies, energyChannel);
            msg.gatherDoubles_(vt,orderChannel);
            ptm.update_(energies);
            msg.scatterDoubles_(ptm.processOrdered_(betas), betaChannel);
            mcenFile << ptm.variableOrdered_(energies);
            mcopFile << ptm.variableOrdered_(vt);
        }



        return 0;
    };
}