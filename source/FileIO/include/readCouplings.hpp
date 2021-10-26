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
    void readNNProfile(base::NNProfile& nnp, std::ifstream& ifp);

    template<typename T, unsigned int latDim, unsigned int spinDim>
    int readLinearInteractions(base::LinearCouplings<T, spinDim>& target, std::string& filename){
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening coupling file " + filename;
                throw std::runtime_error(errormsg); 
            }
            auto string_map = setupCouplingMap();
            auto lat = readLattice<latDim>(filename);
            while (ifp.good()){
                std::string line;
                getline(ifp, line);
                auto entry = readEntryName(line);
                auto it = string_map.find(entry);
                if (it != string_map.end()){
                    switch (string_map.at(entry))
                    {
                    case CouplingType::NN:{
                        auto nnp = base::NNProfile(lat);
                        readNNProfile(nnp,ifp);
                        target.addNN(nnp);
                        break;
                    }
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

    void readNNProfile(base::NNProfile& nnp, std::ifstream& ifp)
    {
        auto read = true;
        while (ifp.good() && read){
            auto string_map = setupNN();
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            auto it = string_map.find(entry);
            if (it != string_map.end()){
                switch (string_map.at(entry))
                {
                case NNParams::BREAK:{
                    read = false;
                    break;
                }
                case NNParams::MAGNITUDE:{
                    nnp.magnitude = readValue<double>(line);
                    break;
                }
                case NNParams::CUTOFF:{
                    nnp.cutoff = readValue<double>(line);
                    break;
                }
                default:
                    break;
                }
            }
        }
    }
}

#endif