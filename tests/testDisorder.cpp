#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/disorderedLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "MonteCarlo/include/VectorModelManager.hpp"
#include "MonteCarlo/include/mcProfile.hpp"
#include <array>
#include <math.h>
#include <iostream>
#include <algorithm>

using namespace classmag;

template<unsigned int dimension>
void printLattice(const geometry::Lattice<dimension> lattice){
    for (unsigned int ii = 0; ii < lattice.n_sites_(); ++ii){
        auto e = lattice.position_(ii);
        for (auto d : e)
            std::cout << d << " ";
        std::cout << "\n";
    }
};

int main(){
    const unsigned int dim = 3;
    auto L = 8u;
    auto runs = 10000u;
    auto occ = 0.75;
    std::array<unsigned int, dim> size;
    size.fill(L);
    auto quotient = 0.0;
    for (auto ii = 0u; ii < runs; ++ii){
        auto sl = geometry::disorderedCubic<dim>(size, occ);
        quotient += static_cast<double>(sl.n_sites_());
    }
    quotient *= 1.0/static_cast<double>(runs);
    quotient *= 1.0/static_cast<double>(pow(L,dim));
    std::cout << quotient/occ << "\n";
    return 0;
};