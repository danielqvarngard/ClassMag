#include <iostream>

#include "Base/include/dipole.hpp"
#include "Geometry/include/predefLattices.hpp"
using namespace classmag;

template<unsigned int dim>
void printEuclidean(geometry::Euclidean<dim> &e){
    for (auto d : e)
        std::cout << d << " ";
    std::cout << "\n";
}

int main(int argc, char* argv[]){
    auto size = std::array<unsigned int,3>{8u, 8u, 8u};
    auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto ep = base::EwaldProfile(lat);
    ep.alpha_ = 1.0;
    ep.magnitude_ = 2.0;
    #if 1
    auto nmin = 2u;
    auto nmax = 30u;
    for (auto ii = nmin; ii <= nmax; ++ii){
        auto trq = 1.0;
        ep.realMirrors_ = ii;
        ep.recMirrors_ = ii;
        ep.alpha_ = 1;
        auto mat  = base::dipoleMatrix(0u, 1u, ep);
        ep.realMirrors_ = ii;
        ep.recMirrors_ = ii;
        trq *= geometry::trace(mat);
        ep.alpha_ = 0.99;
        auto mat2 = base::dipoleMatrix(0u, 1u, ep);
        trq *= 1.0/(geometry::trace(mat2));
        std::cout << 100.0 * abs(1.0 - trq) << "\n"; 
    }
    #endif
}