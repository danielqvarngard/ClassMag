#ifndef CLASSMAG_FILEIO_READLATTICE_HPP
#define CLASSMAG_FILEIO_READLATTICE_HPP

#include <stdexcept>
#include <map>

#include "Geometry/include/predefLattices.hpp"
#include "stringExtractions.hpp"

namespace classmag::fileio{
    enum GeometryOptions{
        PDI_GEO_INVALID = -11,
        PDI_GEO_SIZE = 11,
        PDI_GEO_SUBLATTICE = 12,
        PDI_GEO_BRAVAIS = 13,
        PDI_GEO_DECORATION = 14,
        PDI_GEO_PREDEFINED = 15
    };

    std::map<const std::string, int> mapsetup(){
        std::map<const std::string, int> result;
        result["Size"] = PDI_GEO_INVALID;
        result["Sublattice"] = PDI_GEO_SIZE;
        result["Bravais"] = PDI_GEO_BRAVAIS;
        result["Bravais vectors"] = PDI_GEO_BRAVAIS;
        result["Decoration"] = PDI_GEO_DECORATION;
        result["Predefined"] = PDI_GEO_PREDEFINED;
        result["Predefined lattice"] = PDI_GEO_PREDEFINED;
        return result;
    }

    template<unsigned int latDim>
    geometry::Lattice<latDim> readLattice(const std::string& latticefile){
        auto strmap = mapsetup();
        geometry::Lattice<latDim> lattice;
        try
        {
            std::ifstream ifp("latticefile");
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening lattice file " + latticefile;
                throw std::runtime_error(errormsg); 
            }
            while (ifp.good()){
                std::string str;
                std::stringstream str_stream;
                getline(ifp, str);
                str_stream << str;
                std::string entry;
                auto entry = readEntryName(entry);
                auto it = map.find(str);
                if (it != strmap.end()){
                    switch (strmap.at(str))
                    {
                    case PDI_GEO_SIZE:
                        getline(ifp,str);
                        lattice.setSize_(readRow<unsigned int,latDim>(str));
                        break;
                    case PDI_GEO_BRAVAIS:
                        lattice.setBravais_(readMatrix<double,latDim>(ifp));
                        break;
                    case PDI_GEO_DECORATION:
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
        
        return lattice;
    }
}

#endif