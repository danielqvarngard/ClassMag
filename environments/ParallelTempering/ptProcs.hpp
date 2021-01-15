#ifndef CLASSMAG_ENVIRONMENTS_PTPROCS_HPP
#define CLASSMAG_ENVIRONMENTS_PTPROCS_HPP

#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/StreamManager.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "Parallelism/include/messengers.hpp"
namespace classmag::environments{
    enum MPI_Channels{
        betaChannel = 0, 
        energyChannel = 1, 
        orderChannel = 2,
        invalidChannel = -1
    };

    int mpiPT_hub(unsigned int seed);

    template <unsigned int spinDimension>
    int mpiPT_mc(VectorModelManager<spinDimension> &mc){
        parallelism::Listener msg;
        msg.getDouble_(mc.beta_,betaChannel);
        mc.thermalize_();
        
        for (auto ii = 0u; ii < mc.measurements_(); ++ii){
            mc.update_();
            msg.sendDouble_(mc.energy_(),energyChannel);
            msg.getDouble_(mc.beta_,betaChannel);
        }
        return 0;
    };
}

#endif