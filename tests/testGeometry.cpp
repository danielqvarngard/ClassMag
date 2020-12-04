#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include <array>
#include <iostream>

using namespace classmag::geometry;

int main(){
    auto e1 = Euclidean<3>({0.0,0.0,0.0});
    auto e2 = Euclidean<3>({0.5, 0.5, 0.5});
    auto decoration = std::vector<Euclidean<3>>({e1,e2});
    const std::array<unsigned int,3> systemSize = {4,3,3};
    auto lattice = cubicLattice<3>(systemSize);
    lattice.decorate_(decoration);

    auto n_sites = lattice.n_sites_();

    std::cout << n_sites << "\n";
    

    auto ncs = lattice.neighborCellSites_(0);

    for (auto n : ncs)
        std::cout << n << " ";
    std::cout << "\n";
}