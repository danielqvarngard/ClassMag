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

double volume(const geometry::Lattice<3> &lat){
    auto bravais = lat.getBravais_();
    auto size = lat.getSize_();
    auto V = static_cast<double>(size[0]*size[1]*size[2]);
    V *= abs(bravais[0] * geometry::cross(bravais[1],bravais[2]));
    return V;
}

std::vector<geometry::Euclidean<3>> pyrochlore_ising_axes(){
    auto result = std::vector<geometry::Euclidean<3>>({
        {+1.0, +1.0, +1.0},
        {-1.0, -1.0, +1.0},
        {-1.0, +1.0, -1.0},
        {+1.0, -1.0, -1.0},
        {+1.0, +1.0, +1.0},
        {-1.0, -1.0, +1.0},
        {-1.0, +1.0, -1.0},
        {+1.0, -1.0, -1.0},
        {+1.0, +1.0, +1.0},
        {-1.0, -1.0, +1.0},
        {-1.0, +1.0, -1.0},
        {+1.0, -1.0, -1.0},
        {+1.0, +1.0, +1.0},
        {-1.0, -1.0, +1.0},
        {-1.0, +1.0, -1.0},
        {+1.0, -1.0, -1.0}});
    return result;
}

int main(int argc, char* argv[]){
    auto size = std::array<unsigned int,3>{1u, 1u, 1u};
    auto lat = geometry::pyrochlore(size);
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),volume(lat));
    ep.magnitude_ = 2.0;
    #if 1
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
    for (auto ii = 0u; ii < result.size(); ++ii){
        result[ii] *= 1.0/result[result.size() - 1u];
    }
    print(result);
}