#ifndef CLASSMAG_FILEIO_COUPLINGSTRINGMAPS_HPP
#define CLASSMAG_FILEIO_COUPLINGSTRINGMAPS_HPP

#include <string>
#include <map>
#include "stringExtractions.hpp"

namespace classmag::fileio{
    enum class CouplingType{
        INVALID,
        BREAK,
        NN,
        RKKY,
        DIPOLE
    };

    std::map<const std::string, CouplingType> setupCouplingMap(){
        std::map<const std::string, CouplingType> result;
        result["NN"] = CouplingType::NN;
        result["Nearest neighbor"] = CouplingType::NN;
        result["Nearest neighbour"] = CouplingType::NN;
        result["Nearest Neighbor"] = CouplingType::NN;
        result["Nearest Neighbour"] = CouplingType::NN;
        result["RKKY"] = CouplingType::RKKY;
        result["rkky"] = CouplingType::RKKY;
        result["Dipole"] = CouplingType::DIPOLE;
        return result;
    }

    enum class NNParams{
        INVALID,
        BREAK,
        MAGNITUDE,
        CUTOFF
    };

    std::map<const std::string, NNParams> setupNN(){
        std::map<const std::string, NNParams> result;
        result[breakstring()] = NNParams::BREAK;
        result["Magnitude"] = NNParams::MAGNITUDE;
        result["NN Magnitude"] = NNParams::MAGNITUDE;
        result["NN magnitude"] = NNParams::MAGNITUDE;
        result["Cutoff"] = NNParams::CUTOFF;
        result["Cutoff Radius"] = NNParams::CUTOFF;
        result["Cutoff radius"] = NNParams::CUTOFF;
        return result;
    }

    enum class RKKYParams{
        INVALID,
        BREAK,
        MAGNITUDE,
        KF,
        CUTOFF,
        ALPHA,
        KMIRRORS,
        RMIRRORS
    };

    std::map<const std::string, RKKYParams> setupRKKY(){
        std::map<const std::string,RKKYParams> result;
        result[breakstring()] = RKKYParams::BREAK;
        result["Magnitude"] = RKKYParams::MAGNITUDE;
        result["kf"] = RKKYParams::KF;
        result["Fermi wavevector"] = RKKYParams::KF;
        result["Fermi Wavevector"] = RKKYParams::KF;
        result["Cutoff"] = RKKYParams::CUTOFF;
        result["Cutoff radius"] = RKKYParams::CUTOFF;
        result["Cutoff Radius"] = RKKYParams::CUTOFF;
        result["Alpha"] = RKKYParams::ALPHA;
        result["k Mirrors"] = RKKYParams::KMIRRORS;
        result["k mirrors"] = RKKYParams::KMIRRORS;
        result["Reciprocal Mirrors"] = RKKYParams::KMIRRORS;
        result["Reciprocal mirrors"] = RKKYParams::KMIRRORS;
        result["Real Mirrors"] = RKKYParams::RMIRRORS;
        result["Real mirrors"] = RKKYParams::RMIRRORS; 
        return result;
    }

    enum class DipoleParams{
        INVALID,
        BREAK,
        MAGNITUDE,
        CUTOFF,
        ALPHA,
        KMIRRORS,
        RMIRRORS
    };

    std::map<const std::string, DipoleParams> setupDipole(){
        std::map<const std::string,DipoleParams> result;
        result[breakstring()] = DipoleParams::BREAK;
        result["Magnitude"] = DipoleParams::MAGNITUDE;
        result["Cutoff"] = DipoleParams::CUTOFF;
        result["Cutoff radius"] = DipoleParams::CUTOFF;
        result["Cutoff Radius"] = DipoleParams::CUTOFF;
        result["Alpha"] = DipoleParams::ALPHA;
        result["k Mirrors"] = DipoleParams::KMIRRORS;
        result["k mirrors"] = DipoleParams::KMIRRORS;
        result["Reciprocal Mirrors"] = DipoleParams::KMIRRORS;
        result["Reciprocal mirrors"] = DipoleParams::KMIRRORS;
        result["Real Mirrors"] = DipoleParams::RMIRRORS;
        result["Real mirrors"] = DipoleParams::RMIRRORS; 
        return result;
    }
}

#endif