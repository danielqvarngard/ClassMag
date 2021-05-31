#ifndef CLASSMAG_FILEIO_COUPLINGS_HPP
#define CLASSMAG_FILEIO_COUPLINGS_HPP

#include "Base/include/simulationProcess.hpp"
#include "stringExtractions.hpp"
namespace classmag::fileio{

    enum class CouplingType{
        INVALID,
        NN,
        RKKY,
        DIPOLE
    };

    std::map<const std::string&, CouplingType> setupCouplingMap(){
        std::map<const std::string&, CouplingType> result;
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
        MAGNITUDE,
        CUTOFF
    };

    std::map<const std::string&, NNParams> setupNN(){
        std::map<const std::string&, NNParams> result;
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
        MAGNITUDE,
        KF,
        CUTOFF,
        ALPHA,
        KMIRRORS,
        RMIRRORS
    };

    enum class DipoleParams{
        INVALID,
        MAGNITUDE,
        CUTOFF,
        ALPHA,
        KMIRRORS,
        RMIRRORS
    };

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