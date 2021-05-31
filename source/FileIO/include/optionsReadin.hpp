#ifndef CLASSMAG_FILEIO_OPTIONSREADIN_HPP
#define CLASSMAG_FILEIO_OPTIONSREADIN_HPP

#include <iostream>

#include "Base/include/linearCoupling.hpp"
#include "MonteCarlo/include/HeatBath.hpp"

#include "readLattice.hpp"
#include "readMCO.hpp"
#include "readSpins.hpp"

namespace classmag::fileio{

    enum class OptState{
        INVALID,
        LATTICE,
        LATTICEFILE,
        SPINS,
        SPINFILE,
        MCSETTINGS,
        ASDSETTINGS
    };

    std::map<const std::string, OptState> mapBranch(){
        std::map<const std::string, OptState> result;
        result["Lattice"] = OptState::LATTICE;
        result["Lattice file"] = OptState::LATTICEFILE;
        result["Lattice File"] = OptState::LATTICEFILE;
        result["Spins"] = OptState::SPINS;
        result["Spin file"] = OptState::SPINFILE;
        result["Spin File"] = OptState::SPINFILE;
        result["Monte Carlo settings"] = OptState::MCSETTINGS;
        result["Monte Carlo Settings"] = OptState::MCSETTINGS;
        return result;
    }

    template<typename T, unsigned int dim>
    void printrow(std::array<T, dim> &row){
        for(auto ii = 0u; ii < dim; ++ii)
            std::cout << row[ii] << " ";
    }

    template<typename T, unsigned int dim>
    void printmatrix(std::vector<std::array<T, dim>> &matrix){
        for(auto ii = 0u; ii < matrix.size(); ++ii){
            printrow<T, dim>(matrix[ii]);
            std::cout << "\n";
        }
    }

    template<unsigned int latDim, unsigned int spinDim>
    base::LinearCouplings<spinDim> readModel(const std::string& filename){
        std::ifstream ifp(filename);
        std::string str;
        getline(ifp,str);
        readEntryName(str);

    }

    template<unsigned int latDim, unsigned int spinDim>
    int readOptions(montecarlo::HeatBath<spinDim>& mc, const std::string& filename){

    }
}

#endif