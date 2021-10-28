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

std::vector<std::array<double,8>> comparison_pyrochlore_cubic(){
    auto result = std::vector<std::array<double,8>>({
        {-0.523598775598288, -4.283853877480048, 1.356528775342571, -1.065331984362533, 1.356528775342573, -1.065331984362532, 0.5598676419301706, -0.5235987755982967},
        {-4.283853877480048, -0.523598775598288, -1.065331984362533, 1.356528775342571, -1.065331984362532, 1.356528775342573, -0.5235987755982967, 0.5598676419301706},
        {1.356528775342571, -1.065331984362533, -0.523598775598288, -4.283853877480048, 0.5598676419301706, -0.5235987755982967, 1.356528775342573, -1.065331984362532},
        {-1.065331984362533, 1.356528775342571, -4.283853877480048, -0.523598775598288, -0.5235987755982967, 0.5598676419301706, -1.065331984362532, 1.356528775342573},
        {1.356528775342573, -1.065331984362532, 0.5598676419301706, -0.5235987755982967, -0.523598775598288, -4.283853877480048, 1.356528775342571, -1.065331984362533},
        {-1.065331984362532, 1.356528775342573, -0.5235987755982967, 0.5598676419301706, -4.283853877480048, -0.523598775598288, -1.065331984362533, 1.356528775342571},
        {0.5598676419301706, -0.5235987755982967, 1.356528775342573, -1.065331984362532, 1.356528775342571, -1.065331984362533, -0.523598775598288, -4.283853877480048},
        {-0.5235987755982967, 0.5598676419301706, -1.065331984362532, 1.356528775342573, -1.065331984362533, 1.356528775342571, -4.283853877480048, -0.523598775598288}
    });

    return result; 
}

int main(int argc, char* argv[]){
    auto size = std::array<unsigned int,3>{2u, 2u, 2u};
    auto lat = geometry::pyrochlore(size);
    auto axes = pyrochlore_ising_axes();
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),volume(lat));
    ep.magnitude_ = 0.01;
    auto nmax = 10u;
    ep.realMirrors_ = nmax;
    ep.recMirrors_ = nmax;
    
    for (auto ii = 0u; ii < 16; ++ii){
        for (auto jj = 0u; jj < 16; ++jj){
            auto mat = base::dipoleMatrix(ii, jj, ep);
            auto coupling = axes[ii] * (mat * axes[jj]);
            std::cout << coupling << " ";
        }
        std::cout << "\n";
    }
}