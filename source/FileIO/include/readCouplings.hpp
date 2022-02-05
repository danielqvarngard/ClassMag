#ifndef CLASSMAG_FILEIO_COUPLINGS_HPP
#define CLASSMAG_FILEIO_COUPLINGS_HPP

#include <iostream>

#include "Base/include/simulationProcess.hpp"
#include "Base/include/linearCoupling.hpp"
#include "couplingStringMaps.hpp"
#include "stringExtractions.hpp"
#include "readLattice.hpp"

namespace classmag::fileio{
    void readNNProfile(base::NNProfile& nnp, std::ifstream& ifp);
    void read_dipole_profile(base::DipoleProfile &profile, std::ifstream& ifp);
    void read_rkky_profile(base::RKKYProfile &profile, std::ifstream& ifp);

    template<typename T, unsigned int latDim, unsigned int spinDim>
    int readScalarInteractions(base::LinearCouplings<T, spinDim>& target, const std::string& filename){
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
                    #if 0
                    case CouplingType::DIPOLE:{
                        auto nnp = base::DipoleProfile(lat);
                        read_dipole_profile(nnp,ifp);
                        base::addDipole(target,nnp);
                        break;
                    }
                    #endif
                    case CouplingType::RKKY:{
                        auto nnp = base::RKKYProfile(lat);
                        read_rkky_profile(nnp,ifp);
                        target.addRKKY(nnp);
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

    template<unsigned int latDim>
    int readMatrixInteractions(base::CouplingsMatrixDense& target, const std::string& filename){
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
                    case CouplingType::DIPOLE:{
                        auto nnp = base::DipoleProfile(lat);
                        read_dipole_profile(nnp,ifp);
                        target.addDipole(nnp);
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

    // Can't we cast this into a class of its own and polymorph?
    template<typename T, unsigned int latDim, unsigned int spinDim>
    int readLinearInteractions(
        base::LinearCouplings<T, spinDim>& target, 
        const std::string& filename,
        const geometry::Lattice<latDim>& lat
        ){
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening coupling file " + filename;
                throw std::runtime_error(errormsg); 
            }
            auto string_map = setupCouplingMap();
            //auto lat = readLattice<latDim>(filename);
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
                    case CouplingType::DIPOLE:{
                        auto profile = base::DipoleProfile(lat);
                        read_dipole_profile(profile, ifp);
                        target.addDipole(profile);
                        break;
                    }
                    case CouplingType::RKKY:{
                        auto profile = base::RKKYProfile(lat);
                        read_rkky_profile(profile, ifp);
                        target.addRKKY(profile);
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

    void read_dipole_profile(base::DipoleProfile& profile, std::ifstream& ifp)
    {
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
                    profile.magnitude_ = readValue<double>(line);
                    break;
                }
                case DipoleParams::ALPHA:{
                    profile.alpha_ = readValue<double>(line);
                    alpha_set = true;
                    break;
                }
                case DipoleParams::KMIRRORS:{
                    profile.recMirrors_ = readValue<double>(line);
                    break;
                }
                case DipoleParams::RMIRRORS:{
                    profile.realMirrors_ = readValue<double>(line);
                    break;
                }
                default:
                    break;
                }
            }

            if (!alpha_set){
                profile.alpha_ = base::optimAlpha(
                    profile.lattice_.n_sites_(), 
                    base::volume(profile.lattice_));
            }
        }
    }

    void read_rkky_profile(base::RKKYProfile &profile, std::ifstream& ifp){
        auto read = true;
        while (ifp.good() && read){
            auto string_map = setupRKKY();
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            auto it = string_map.find(entry);
            bool alpha_set = false;
            if (it != string_map.end()){
                switch (string_map.at(entry))
                {
                case RKKYParams::BREAK:{
                    read = false;
                    break;
                }
                case RKKYParams::MAGNITUDE:{
                    profile.magnitude = readValue<double>(line);
                    break;
                }
                case RKKYParams::KF:{
                    profile.k_F = readValue<double>(line);
                    break;
                }
                case RKKYParams::CUTOFF:{
                    profile.cutoff = readValue<double>(line);
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