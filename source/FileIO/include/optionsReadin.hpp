#ifndef CLASSMAG_FILEIO_OPTIONSREADIN_HPP
#define CLASSMAG_FILEIO_OPTIONSREADIN_HPP

#include <iostream>
#include <fstream>
#include <map>
//#include "Base/include/linearCoupling.hpp"
//#include "MonteCarlo/include/HeatBath.hpp"

//#include "readLattice.hpp"
#include "readMCO.hpp"
//#include "readSpins.hpp"
#include "scancmd.hpp"

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

    inline std::map<const std::string, OptState> mapBranch(){
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

    std::string get_input_file_name(int argc, char* argv[]);

    std::ifstream get_input_stream(int argc, char* argv[]);

}

#endif