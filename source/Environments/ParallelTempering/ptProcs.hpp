#ifndef CLASSMAG_ENVIRONMENTS_PTPROCS_HPP
#define CLASSMAG_ENVIRONMENTS_PTPROCS_HPP

#include <iostream>
#include <mpi.h>

#include <utility>

#include "PT_Return.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/StreamManager.hpp"
#include "MonteCarlo/include/HeatBath.hpp"
#include "MonteCarlo/include/XYModelManager.hpp"
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

    struct PT_Input{
        int measurement_runs_;
        int thermalization_runs_;
        std::vector<double> betas_;
    };

    void MPI_PT_Hub_Thermalize(const PT_Input& input);

    void MPI_PT_Hub(const std::vector<double>& betas, const PT_Input& input);

    PT_Return MPI_PT_Hub(const PT_Input& input);

    PT_Return MPI_PT_Hub(const PT_Input& input, const unsigned int orderCount);

    template <unsigned int spinDimension> // Eww, we should use polymorphism or something
    int mpiPT_mc(
        montecarlo::XYModelManager<spinDimension> &mc,
        montecarlo::VectorModel_Profile& mco){
        parallelism::EdgeNodeMessenger msg;

        for (auto ii = 0; ii < mco.thermalization_; ++ii){
            msg.getDouble_(mc.beta_, betaChannel);
            mc.update_();
        }
        

        for (auto ii = 0u; ii < mco.measurement_; ++ii){
            mc.update_();
            auto e = mc.energy_();
            msg.sendDouble_(e, energyChannel);
            msg.getDouble_(mc.beta_, betaChannel);
        }
        return 0;
    };

    template <unsigned int spinDimension> // Eww, we should use polymorphism or something
    int mpiPT_mc(
        montecarlo::HeatBath<spinDimension> &mc,
        montecarlo::VectorModel_Profile& mco){
        parallelism::EdgeNodeMessenger msg;
        msg.getDouble_(mc.beta_, betaChannel);
        mc.update_(mco.thermalization_);
        for (auto ii = 0u; ii < mco.measurement_; ++ii){
            mc.update_();
            auto e = mc.energy_();
            msg.sendDouble_(e, energyChannel);
            msg.getDouble_(mc.beta_, betaChannel);
        }
        return 0;
    };

    template <unsigned int spinDimension>
    int mpiPT_mc(montecarlo::ClockModelManager<spinDimension> &mc){
        parallelism::EdgeNodeMessenger msg;
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