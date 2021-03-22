#include "include/hasxStagMag.hpp"

namespace classmag::base{
    std::pair<OrderParameter<3>, OrderParameter<3>> chasxStagMag_split(
    const geometry::Lattice<3> &lattice){
        
        auto partitions = lattice.partitions_();
        std::function<double(const SpinStructure<3>&)> 
            shellOrder = 
            [partitions](const SpinStructure<3> &spins){
            auto m = geometry::Euclidean<3>({0.0, 0.0, 0.0});
            for (auto b = 0u; b < 2u; ++b){
                auto staggerFactor = 2.0 * static_cast<double>(b % 2) - 1.0;
                for (unsigned int site = partitions[b]; site < partitions[b + 1]; ++site){
                    m += staggerFactor * spins[site];
                }
            }
            return geometry::norm(m);
        };

        std::function<double(const SpinStructure<3>&)> 
            centerOrder = 
            [partitions](const SpinStructure<3> &spins){
            auto m = geometry::Euclidean<3>({0.0, 0.0, 0.0});
            for (auto b = 2u; b < partitions.size() - 1; ++b){
                auto staggerFactor = 2.0 * static_cast<double>(b % 2) - 1.0;
                for (unsigned int site = partitions[b]; site < partitions[b + 1]; ++site){
                    m += staggerFactor * spins[site];
                }
            }
            return geometry::norm(m);
        };

        auto opShell = OrderParameter<3>("Shell order", shellOrder);
        auto opCenter = OrderParameter<3>("Center order", centerOrder);
        auto opPair = std::pair<OrderParameter<3>,OrderParameter<3>>({opShell,opCenter});
        return opPair;
    };

    std::pair<OrderParameter<3>, OrderParameter<3>> chasxMag_split(
    const geometry::Lattice<3> &lattice){
        
        auto partitions = lattice.partitions_();
        std::function<double(const SpinStructure<3>&)> 
            shellOrder = 
            [partitions](const SpinStructure<3> &spins){
            auto m = geometry::Euclidean<3>({0.0, 0.0, 0.0});
                for (unsigned int site = partitions[0]; site < partitions[2]; ++site){
                    m += spins[site];
                }
            return geometry::norm(m);
        };

        std::function<double(const SpinStructure<3>&)> 
            centerOrder = 
            [partitions](const SpinStructure<3> &spins){
            auto m = geometry::Euclidean<3>({0.0, 0.0, 0.0});
                for (unsigned int site = partitions[2]; site < partitions[3]; ++site){
                    m += spins[site];
                }
            return geometry::norm(m);
        };

        auto opShell = OrderParameter<3>("Shell order", shellOrder);
        auto opCenter = OrderParameter<3>("Center order", centerOrder);
        auto opPair = std::pair<OrderParameter<3>,OrderParameter<3>>({opShell,opCenter});
        return opPair;
    };
}