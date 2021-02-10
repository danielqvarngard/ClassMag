#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/predefLattices.hpp"
#include <array>
#include <iostream>

using namespace classmag::geometry;

template<unsigned int dimension>
void printLattice(const Lattice<dimension> lattice){
    for (unsigned int ii = 0; ii < lattice.n_sites_(); ++ii){
        auto e = lattice.position_(ii);
        for (auto d : e)
            std::cout << d << " ";
        std::cout << "\n";
    }
};

int main(){

    const auto systemSize = std::array<unsigned int,3>({2,2,2});
    auto lattice = has100(systemSize);
    
    printLattice(lattice);
};