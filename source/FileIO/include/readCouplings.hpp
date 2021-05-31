#ifndef CLASSMAG_FILEIO_COUPLINGS_HPP
#define CLASSMAG_FILEIO_COUPLINGS_HPP

#include "Base/include/simulationProcess.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/rkky.hpp"
#include "Base/include/dipole.hpp"
#include "couplingStringMaps.hpp"
#include "stringExtractions.hpp"
#include "readLattice.hpp"

namespace classmag::fileio{

    template<unsigned int latDim>
    void readNNProfile(base::NNProfile<latDim>& nnp, std::ifstream& ifp){

    }

    template<unsigned int latDim, unsigned int spinDim>
    int readLinearInteractions(base::LinearCouplings<spinDim>& target, std::string& filename){
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening coupling file " + filename;
                throw std::runtime_error(errormsg); 
            }
            auto strmap = setupCouplingMap();
            auto lat = readLattice<latDim>(filename);
            while (ifp.good()){
                std::string str;
                std::stringstream str_stream;
                getline(ifp, str);
                str_stream << str;
                std::string entry;
                auto entry = readEntryName(entry);
                auto it = strmap.find(str);
                if (it != strmap.end()){
                    switch (strmap.at(str))
                    {
                    case CouplingType::NN:
                        auto nnp = base::NNProfile(lat);
                        readNNProfile(nnp,ifp);
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            abort();
        }
        
        return 0;
    };
}

#endif