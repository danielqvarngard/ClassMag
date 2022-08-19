#ifndef CLASSMAG_FILEIO_READDEFAULTSETTINGS_HPP
#define CLASSMAG_FILEIO_READDEFAULTSETTINGS_HPP

#include <iostream>
#include <fstream>

#include "couplingStringMaps.hpp"
#include "filesystem.hpp"
#include "stringExtractions.hpp"
#include "readLattice.hpp"

namespace classmag::fileio{
    void readNNProfile(DefaultMonteCarloParameters& x, std::ifstream& ifp);
    void read_dipole_profile(DefaultMonteCarloParameters& x, std::ifstream& ifp);
    void read_rkky_profile(DefaultMonteCarloParameters& x, std::ifstream& ifp);

    enum class MCDefaultOptions{
        INVALID,
        BREAK,
        NN,
        RKKY,
        DIPOLE,
        THERMALIZATIONS,
        SKIPS,
        MEASUREMENTS,
        OVERRELAXATIONS,
        RESAMPLES,
        ORDERPARAMETERS,
        SIZE,
        TEMPERATURES
    };

    std::map<const std::string, MCDefaultOptions> setupHackyMap(){
        std::map<const std::string, MCDefaultOptions> result;
        result["NN"] = MCDefaultOptions::NN;
        result["Nearest neighbor"] = MCDefaultOptions::NN;
        result["Nearest neighbour"] = MCDefaultOptions::NN;
        result["Nearest Neighbor"] = MCDefaultOptions::NN;
        result["Nearest Neighbour"] = MCDefaultOptions::NN;
        result["RKKY"] = MCDefaultOptions::RKKY;
        result["rkky"] = MCDefaultOptions::RKKY;
        result["Dipole"] = MCDefaultOptions::DIPOLE;
        result["dipole"] = MCDefaultOptions::DIPOLE;
        result["Temperatures"] = MCDefaultOptions::TEMPERATURES;
        result["Thermalizations"] = MCDefaultOptions::THERMALIZATIONS;
        result["Skips"] = MCDefaultOptions::SKIPS;
        result["Measurements"] = MCDefaultOptions::MEASUREMENTS;
        result["Overrelaxations"] = MCDefaultOptions::OVERRELAXATIONS;
        result["Resamples"] = MCDefaultOptions::RESAMPLES;
        result["Order parameters"] = MCDefaultOptions::ORDERPARAMETERS;
        result["Size"] = MCDefaultOptions::SIZE;
        result["System size"] = MCDefaultOptions::SIZE;
        return result;
    };


    DefaultMonteCarloParameters ReadDefaultMonteCarloOptions
    (
        const std::string& filename
    ){
        DefaultMonteCarloParameters result;
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening coupling file " + filename;
                throw std::runtime_error(errormsg); 
            }
            auto string_map = setupHackyMap();
            while (ifp.good()){
                std::string line;
                getline(ifp, line);
                auto entry = readEntryName(line);
                auto it = string_map.find(entry);
                if (it != string_map.end()){
                    switch (string_map.at(entry))
                    {
                    case MCDefaultOptions::NN:{
                        readNNProfile(result,ifp);
                        break;
                    }
                    case MCDefaultOptions::DIPOLE:{
                        read_dipole_profile(result, ifp);
                        break;
                    }
                    case MCDefaultOptions::RKKY:{
                        read_rkky_profile(result, ifp);
                        break;
                    }
                    case MCDefaultOptions::THERMALIZATIONS:{
                        result.thermalizations = readValue<int>(line);
                        break;
                    }
                    case MCDefaultOptions::MEASUREMENTS:{
                        result.measurements = readValue<int>(line);
                        break;
                    }
                    case MCDefaultOptions::SKIPS:{
                        result.skipped_sweeps = readValue<int>(line);
                        break;
                    }
                    case MCDefaultOptions::SIZE:{
                        result.system_size = read_list<int>(line);
                        break;
                    }
                    case MCDefaultOptions::TEMPERATURES:{
                        result.temperatures = read_list<double>(line);
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
        
        return result;
    };

    void readNNProfile(DefaultMonteCarloParameters& x, std::ifstream& ifp)
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
                    x.nearest_neighbor_strength = readValue<double>(line);
                    break;
                }
                case NNParams::CUTOFF:{
                    x.nearest_neighbor_cutoff = readValue<double>(line);
                    break;
                }
                default:
                    break;
                }
            }
        }
    };

    void read_dipole_profile(DefaultMonteCarloParameters& x, std::ifstream& ifp){
        auto read = true;
        while (ifp.good() && read){
            auto string_map = setupDipole();
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            auto it = string_map.find(entry);
            bool alpha_set = false;
            if (it != string_map.end()){
                switch (string_map.at(entry))
                {
                case DipoleParams::BREAK:{
                    read = false;
                    break;
                }
                case DipoleParams::MAGNITUDE:{
                    x.dipole_strength = readValue<double>(line);
                    break;
                }
                case DipoleParams::ALPHA:{
                    x.dipole_alpha = readValue<double>(line);
                    alpha_set = true;
                    break;
                }
                case DipoleParams::KMIRRORS:{
                    x.dipole_k_reciprocal = readValue<double>(line);
                    break;
                }
                case DipoleParams::RMIRRORS:{
                    x.dipole_k_real = readValue<double>(line);
                    break;
                }
                default:
                    break;
                }
            }
        }
    }

    void read_rkky_profile(DefaultMonteCarloParameters& x, std::ifstream& ifp){
        auto read = true;
        while (ifp.good() && read){
            auto string_map = setupRKKY();
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            auto it = string_map.find(entry);
            if (it != string_map.end()){
                switch (string_map.at(entry))
                {
                case RKKYParams::BREAK:{
                    read = false;
                    break;
                }
                case RKKYParams::MAGNITUDE:{
                    x.rkky_strength = readValue<double>(line);
                    break;
                }
                case RKKYParams::KF:{
                    x.rkky_wavevector = readValue<double>(line);
                    break;
                }
                case RKKYParams::CUTOFF:{
                    x.rkky_cutoff = readValue<double>(line);
                    break;
                }
                default:
                    break;
                }
            }
        }
    };
}

#endif