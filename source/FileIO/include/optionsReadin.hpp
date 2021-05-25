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

    enum GeometryOptions{
        PDI_GEO_INVALID = -11,
        PDI_GEO_SIZE = 11,
        PDI_GEO_SUBLATTICE = 12,
        PDI_GEO_BRAVAIS = 13,
        PDI_GEO_DECORATION = 14,
        PDI_GEO_PREDEFINED = 15
    };

    enum MontecarloOptions{
        PDI_MCO_INVALID = -21,
        PDI_MCO_THERMRUNS = 21,
        PDI_MCO_MEASRUNS = 22,
        PDI_MCO_SKIPS = 23
    };

    std::string readEntryName(const std::string &row){
        std::stringstream stream;
        stream << row;
        bool finito = false;
        std::string entry;
        while (!stream.eof() && !finito){
            
            std::string temp;
            stream >> temp;
            if (temp.compare("=") != 0){
                if (!entry.empty())
                    entry += " ";
                entry += temp;
            }
            else
                finito = true;
        }
        return entry;
    }

    int readDelimiter(const std::string &line){
        if (line.find('{') < line.npos)
            return -1;
        if (line.find('}') < line.npos)
            return +1;
        return 0;
    }

    bool matrixEnd(const std::string &line){
        if (line.find('}') < line.npos)
            return true;
        else
            return false;
    }

    template<typename T, unsigned int dim>
    std::array<T, dim> readRow(const std::string line){
        std::stringstream stream;
        std::array<T, dim> result;
        stream << line;
        unsigned int ii = 0;
        while (!stream.eof()){
            std::string temp_str;
            T value;
            stream >> temp_str;
            if (std::stringstream(temp_str) >> value){
                result[ii] = value;
                ++ii;
            }
        }
        if (ii < dim)
            std::cout << "error reading entries\n";
        return result;
        
    }

    template<typename T, unsigned int columns>
    std::vector<std::array<T, columns>> readMatrix(std::ifstream &ifp){
        std::vector<std::array<T,columns>> result;
        bool endOfMatrix = false;
        while (!endOfMatrix){
            std::string line;
            getline(ifp, line);
            if (matrixEnd(line) || ifp.eof())
                endOfMatrix = true;
            else{
                auto row = readRow<T,columns>(line);
                result.push_back(row);
            }
        }
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
}

#endif