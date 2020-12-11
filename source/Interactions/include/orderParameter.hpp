#ifndef CLASSMAG_MODELS_ORDERPARAMETER_HPP
#define CLASSMAG_MODELS_ORDERPARAMETER_HPP

#include <functional>
#include <utility>
#include <numeric>
#include <string>

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/lattice.hpp"
#include "spinStructure.hpp"

namespace classmag::models{
    template <unsigned int spinDimension>
    class OrderParameter{
        public:
        OrderParameter(const std::string &name, 
            std::function <double(const SpinStructure<spinDimension>&)> 
                orderParameter):
            name_(name),
            measure_(orderParameter){

        }

        const std::string name_;

        std::function
            <double(const SpinStructure<spinDimension>&)>
            measure_;
    };

    template <unsigned int spinDimension>
    class OrderParameters : public std::vector<OrderParameter<spinDimension>>{
        
    };

    template <unsigned int spinDimension>
    OrderParameter<spinDimension> magnetization(){
        auto measurement = [](const SpinStructure<spinDimension> &spins){
            auto sum = spins[0];
            for (unsigned int ii = 1; ii < spins.size(); ++ii)
                sum += spins[ii];

            auto magnitude = geometry::norm(sum);
            return magnitude;
        };
        auto name = "Magnetization";
        auto op = OrderParameter<3>(name, measurement);
        return op;
    };
    
    template<unsigned int dim, unsigned int spinDimension>
    OrderParameter<spinDimension> clusterOrder(
        const geometry::Lattice<dim> lattice){
            std::function<double(const SpinStructure<spinDimension>& spins)> order =
                [lattice](const SpinStructure<spinDimension>& spins){
                auto result = 0.0;
                for (unsigned int site = 0; site < spins.size(); ++site){
                    auto loopIndices = lattice.neighborCellSites_(site);
                    auto n_neighbors = static_cast<double>(loopIndices.size());
                    for (auto neighborSite : loopIndices)
                        result += spins[site]*spins[neighborSite]/n_neighbors;
                }
                auto n_sites = static_cast<double>(spins.size());
                result /= n_sites;
                return result;
            };

            auto clusterFunction = OrderParameter<spinDimension>("Cluster",order);
            return clusterFunction;
        };

    std::pair<geometry::Lattice<3>, OrderParameter<3>> DefaultTsaiSim(
        const std::array<unsigned int,3> &systemSize);
}

#endif