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
                std::string line;
                getline(ifp, line);
                auto entry = readEntryName(line);
                auto it = strmap.find(entry);
                if (it != strmap.end()){
                    switch (strmap.at(entry))
                    {
                    case CouplingType::NN:
                        auto nnp = base::NNProfile(lat);
                        readNNProfile(nnp,ifp);
                        target.addNN(nnp);
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

    void readNNProfile(base::NNProfile& nnp, std::ifstream& ifp)
    {
        auto read = true;
        while (ifp.good() && read){
            auto strmap = setupNN();
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            auto it = strmap.find(entry);
            if (it != strmap.end()){
                switch (strmap.at(entry))
                {
                case NNParams::BREAK:
                    read = false;
                    break;
                case NNParams::MAGNITUDE:
                    nnp.magnitude = readValue<double>(line);
                    break;
                case NNParams::CUTOFF:
                    nnp.cutoff = readValue<double>(line);
                    break;
                default:
                    break;
                }
            }
        }
    }
}

#endif