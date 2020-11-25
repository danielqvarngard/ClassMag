#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/lattice.hpp"
#include <array>
#include <iostream>

using namespace classmag::geometry;

int main(){
    auto e1 = Euclidean<2>({1.0, 2.0});
    auto e2 = Euclidean<2>({0.0, 3.0});
    std::cout << norm(proj(e1,e2)) << "\n";
}