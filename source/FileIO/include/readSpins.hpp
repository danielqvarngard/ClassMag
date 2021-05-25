#ifndef CLASSMAG_FILEIO_READSPINS_HPP
#define CLASSMAG_FILEIO_READSPINS_HPP

#include "Base/include/simulationProcess.hpp"
#include "stringExtractions.hpp"
namespace classmag::fileio{
    template<unsigned int spinDim>
    int readSpins(
        std::ifstream& spinfile, 
        base::SimulationBase<spinDim>& sim)
    {
        
        base::SpinStructure<spindim> structure;
        sim.setSpinstructure_();
    }
}

#endif