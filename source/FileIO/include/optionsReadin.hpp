#ifndef CLASSMAG_FILEIO_OPTIONSREADIN_HPP
#define CLASSMAG_FILEIO_OPTIONSREADIN_HPP

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>

#include "Base/include/linearCoupling.hpp"

#include "readLattice.hpp"
#include "readMCO.hpp"
#include "readSpins.hpp"

namespace classmag::fileio{

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
}

#endif