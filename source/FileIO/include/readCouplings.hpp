#ifndef CLASSMAG_FILEIO_COUPLINGS_HPP
#define CLASSMAG_FILEIO_COUPLINGS_HPP

#include "Base/include/simulationProcess.hpp"
#include "couplingStringMaps.hpp"
#include "stringExtractions.hpp"

namespace classmag::fileio{

    template<unsigned int latDim, unsigned int spinDim>
    int readLinearInteractions(base::LinearCouplings<spinDim>& target, std::string& filename){
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening coupling file " + filename;
                throw std::runtime_error(errormsg); 
            }
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
                    case GeoState::SIZE:
                        getline(ifp,str);
                        lattice.setSize_(readRow<unsigned int,latDim>(str));
                        break;
                    case GeoState::BRAVAIS:
                        lattice.setBravais_(readMatrix<double,latDim>(ifp));
                        break;
                    case GeoState::DECORATION:
                        lattice.decorate_(readMatrix<double,latDim>(ifp));
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