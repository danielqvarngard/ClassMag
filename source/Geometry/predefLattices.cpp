#include "include/predefLattices.hpp"

namespace classmag::geometry{
    Lattice<2> squareLattice(const std::array<unsigned int,2> &systemSize){
        auto e1 = Euclidean<2>({1.0,0.0});
        auto e2 = Euclidean<2>({0.0,1.0});
        const auto bravais = std::array<Euclidean<2>,2> ({e1, e2});
        auto result = Lattice<2>(bravais,systemSize);
        return result;
    };
}