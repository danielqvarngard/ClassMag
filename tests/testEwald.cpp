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
    auto size = std::array<unsigned int,3>{2u, 2u, 2u};
    auto lat = geometry::cubicLattice<3>(size);
    auto ep = base::EwaldProfile(lat);
    ep.alpha_ = 1.0;
    ep.magnitude_ = 2.0;
    #if 1
    auto nmin = 1u;
    auto nmax = 5u;
    for (auto n = nmin; n < nmax; ++n){
        ep.realMirrors_ = n;
        ep.recMirrors_ = n;
        auto mat = base::dipoleMatrix(0u,1u,ep);
        std::cout << mat[0][0] << " ";
    }
    std::cout << "\n";
    #endif
}