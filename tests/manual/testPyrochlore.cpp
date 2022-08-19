#include <iostream>

#include "Base/include/dipole.hpp"
#include "Base/include/dipole_spherical.hpp"
#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "Geometry/include/predefLattices.hpp"
using namespace classmag;

template<unsigned int dim>
void print(geometry::Euclidean<dim> &e){
    for (auto d : e)
        std::cout << d << " ";
    std::cout << "\n";
}

void print(geometry::Matrix<3,3>& m)
{
    for (auto e : m)
        print<3>(e);
}

void print(std::vector<double> &v){
    for (auto d : v)
        std::cout << d << " ";
    std::cout << "\n";
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
    auto lat = geometry::Lattice<3>(geometry::cubicLattice<3>(size));
    auto axis = geometry::Euclidean<3>({1.0, 0.0, 0.0});
    auto ep = base::DipoleProfile(lat);
    ep.length_ = 2.0;
    auto pata_couplings = comparison_pyrochlore_cubic();
    //ep.alpha_ = base::optimAlpha(lat.n_sites_(),volume(lat));
    ep.alpha_ = 2.0;
    ep.magnitude_ = 1.0;
    auto nmax = 20u;
    ep.realMirrors_ = nmax;
    ep.recMirrors_ = nmax;

    #if 1
    for (auto ii = 0u; ii < 8u; ++ii){
        for (auto jj = 0u; jj < 8u; ++jj){
            auto mat = base::dipole_spherical_matrix(ii, jj, ep);
            //mat = mat + (4.0 * base::pi()) * geometry::eye<3>();
            auto coupling = axis * (mat * axis);
            std::cout << coupling/pata_couplings[ii][jj] << " ";
        }
        std::cout << "\n";
    }
    #endif
}