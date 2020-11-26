#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/lattice.hpp"
#include <array>
#include <iostream>

using namespace classmag::geometry;

int main(){
    auto e1 = Euclidean<2>({1.0, 0.0});
    auto e2 = Euclidean<2>({1.0, 1.0});
    e2 -= proj(e2,e1);

    auto bravais = std::array<Euclidean<2>,2>({e1,e2});
    auto systemSize = std::array<unsigned int,2>({3,3});
    Lattice<2> lattice(bravais,systemSize);

    auto n_sites = systemSize[0]*systemSize[1];
    for (unsigned int ii = 0; ii < n_sites; ++ii)
        std::cout << lattice.position_(ii)[0] << " " << lattice.position_(ii)[1] << "\n";
    
}