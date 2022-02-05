#ifndef CLASSMAG_FILEIO_READLATTICE_HPP
#define CLASSMAG_FILEIO_READLATTICE_HPP

#include <stdexcept>
#include <map>

#include "Geometry/include/predefLattices.hpp"
#include "stringExtractions.hpp"

namespace classmag::fileio{
    enum class GeoState{
        INVALID,
        SIZE,
        SUBLATTICE,
        BRAVAIS,
        DECORATION,
        PREDEFINED
    };

    inline std::map<const std::string, GeoState> mapsetup(){
        std::map<const std::string, GeoState> result;
        result["Size"] = GeoState::SIZE;
        result["Sublattice"] = GeoState::SUBLATTICE;
        result["Bravais"] = GeoState::BRAVAIS;
        result["Bravais vectors"] = GeoState::BRAVAIS;
        result["Decoration"] = GeoState::DECORATION;
        result["Predefined"] = GeoState::PREDEFINED;
        result["Predefined lattice"] = GeoState::PREDEFINED;
        return result;
    }

    template<unsigned int latDim>
    geometry::Lattice<latDim> readLattice(const std::string& latticefile){
        auto strmap = mapsetup();
        geometry::Lattice<latDim> lattice;
        try
        {
            std::ifstream ifp(latticefile);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening lattice file " + latticefile;
                throw std::runtime_error(errormsg); 
            }
            while (ifp.good()){
                std::string str;
                getline(ifp, str);
                auto entry = readEntryName(str);
                
                auto it = strmap.find(entry);
                if (it != strmap.end()){
                    switch (strmap.at(entry))
                    {
                    case GeoState::SIZE:{
                        lattice.setSize_(readRow<unsigned int,latDim>(str));
                        break;
                    }
                    case GeoState::BRAVAIS:{
                        lattice.setBravais_(readMatrix<double,latDim>(ifp));
                        break;
                    }
                    case GeoState::DECORATION:{
                        lattice.decorate_(readMatrix<double,latDim>(ifp));
                        break;
                    }
                    case GeoState::PREDEFINED:{
                        
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
        
        return lattice;
    }

    template<unsigned int latDim>
    geometry::Lattice<latDim> readPredefLattice(){

    }

    template<unsigned int latDim>
    std::array<unsigned int, latDim> readSize(const std::string& latticefile){
        auto strmap = mapsetup();
        std::array<unsigned int, latDim> size;
        try
        {
            std::ifstream ifp(latticefile);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening lattice file " + latticefile;
                throw std::runtime_error(errormsg); 
            }
            while (ifp.good()){
                std::string str;
                getline(ifp, str);
                auto entry = readEntryName(str);
                auto it = strmap.find(entry);
                if (it != strmap.end()){
                    switch (strmap.at(entry))
                    {
                    case GeoState::SIZE:{
                        size = readRow<unsigned int,latDim>(str);
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
        
        return size;
    }
}

#endif