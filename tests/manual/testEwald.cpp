#include <iostream>

#include "Base/include/dipole.hpp"
#include "Geometry/include/predefLattices.hpp"
using namespace classmag;

template<unsigned int dim>
void print(geometry::Euclidean<dim> &e){
    for (auto d : e)
        std::cout << d << " ";
    std::cout << "\n";
}

void print(std::vector<double> &v){
    for (auto d : v)
        std::cout << d << " ";
    std::cout << "\n";
}

template<unsigned int m, unsigned int n>
void print(geometry::Matrix<m,n>& mat){
    for (auto v : mat)
        print(v);
    
}

int main(int argc, char* argv[]){
    auto size = std::array<unsigned int,3>{8u, 3u, 8u};
    auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),base::volume(lat));
    ep.magnitude_ = 2.0;
    #if 0
    auto nmin = 2u;
    auto nmax = 10u;
    std::vector<double> result;
    for (auto ii = nmin; ii <= nmax; ++ii){
        ep.realMirrors_ = ii;
        ep.recMirrors_ = ii;
        auto mat = base::dipoleMatrix(0u, 0u, ep);
        auto trq = geometry::trace(mat);
        result.push_back(trq);
    }
    #endif
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    auto mat = base::dipoleMatrix(0u,0u,ep);
    print(mat);

}