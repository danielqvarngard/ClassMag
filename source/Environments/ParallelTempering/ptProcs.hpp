#ifndef CLASSMAG_ENVIRONMENTS_PTPROCS_HPP
#define CLASSMAG_ENVIRONMENTS_PTPROCS_HPP

#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/StreamManager.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/ClockModelManager.hpp"
#include "MonteCarlo/include/extendedEnsembles.hpp"
#include "Parallelism/include/messengers.hpp"

namespace classmag::environments{
    enum MPI_Channels{
        betaChannel = 0, 
        energyChannel = 1, 
        orderChannel = 2,
        invalidChannel = -1
    };

    int mpiPT_hub(const std::vector<double> &temperatures, const unsigned int measurementCount);

    int mpiPT_hub(
        const std::vector<double> &temperatures, 
        const unsigned int measurementCount,
        const unsigned int orderCount);

    template <unsigned int spinDimension>
    int mpiPT_mc(
        montecarlo::VectorModelManager<spinDimension> &mc){
        parallelism::Listener msg;
        msg.getDouble_(mc.beta_, betaChannel);
        mc.thermalize_();
        
        for (auto ii = 0u; ii < mc.measurements_(); ++ii){
            mc.update_();
            msg.sendDouble_(mc.energy_(), energyChannel);
            msg.sendDouble_(mc.measure_(), orderChannel);
            msg.getDouble_(mc.beta_, betaChannel);
        }
        return 0;
    };

    template <unsigned int spinDimension>
    int mpiPT_mc(montecarlo::ClockModelManager<spinDimension> &mc){
        parallelism::Listener msg;
        msg.getDouble_(mc.beta_, betaChannel);
        mc.thermalize_();
        
        for (auto ii = 0u; ii < mc.measurements_(); ++ii){
            mc.update_();
            msg.sendDouble_(mc.energy_(), energyChannel);
            msg.sendDouble_(mc.measure_(), orderChannel);
            msg.getDouble_(mc.beta_, betaChannel);
        }
        return 0;
    };
}

#endif