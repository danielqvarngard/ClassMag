#include "include/default_orderparameters.hpp"

namespace classmag::base
{
    OrderParameter<3> MagnetizationElement(unsigned int cartesian_index){
        std::function<double(const SpinStructure<3>&)> element =
            [cartesian_index](const SpinStructure<3>& spins){
                auto result = 0.0;
                for (auto s : spins)
                    result += s[cartesian_index];
                return result;
            };
        
        std::string name = "M_";
        name += std::to_string(cartesian_index);
        auto result = OrderParameter(name, element);
        return result;
    };

    std::vector<OrderParameter<3>> Magnetization()
    {
        std::vector<OrderParameter<3>> result;
        for (auto ii = 0; ii < 3; ++ii){
            result.push_back(MagnetizationElement(ii));
        }
        return result;
    }

    OrderParameter<3> BipartiteStaggeredMagnetization(
        const std::vector<unsigned int> sublattice_partition_indices
    ){
        std::function<double(const SpinStructure<3>&)> staggered_magnetization = 
            [sublattice_partition_indices](const SpinStructure<3>& spins){

            auto m = geometry::Euclidean<3>({0.0, 0.0, 0.0});
            for (auto b = 0u; b < 2u; ++b)
            {
                auto stagger_factor = 2.0 * static_cast<double>(b % 2) - 1.0;
                for (unsigned int site = 
                    sublattice_partition_indices[b];
                    site < sublattice_partition_indices[b + 1]; ++site)
                {
                    m += stagger_factor * spins[site];
                }
            }
            return geometry::norm(m);
        };

        std::string name = "Staggered magnetization";
        auto result = OrderParameter<3>(name, staggered_magnetization);
        return result;
    }
}