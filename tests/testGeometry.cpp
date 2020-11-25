#include "Geometry/include/qartesian.hpp"
#include "Geometry/include/lattice.hpp"
#include <array>
#include <iostream>

using namespace classmag::geometry;

int main(){
    const Qartesian<2> a0 = {1.0, 0.0};
    const Qartesian<2> a1 = {0.0, 1.0};
    const std::array<Qartesian<2>,2> bravais = {a0, a1};
    const std::array<unsigned int,2> systemSize = {3, 3};
    Lattice<2> lattice(bravais,systemSize);

    for (unsigned int ii = 0; ii < 9; ++ii)
    {
        auto q = lattice.position_(ii);
        std::cout << q.get(0) << " " << q.get(1) << "\n";
    }
}