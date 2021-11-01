#include <iostream>
#include <sstream>
#include <vector>


#include "Geometry/include/predefLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/dipole.hpp"
#include "Base/include/linearCoupling.hpp"
#include "MonteCarlo/include/HeatBath.hpp"
#include "MonteCarlo/include/ClockModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;

class mock{
    public:
    unsigned int n_sites = 6u;
    
    inline unsigned int columnIndex(const unsigned int a) const
    {
        return (a * n_sites - ((a + 1u) * a)/2u);
    }
    inline unsigned int index(const unsigned int a, const unsigned int b) const 
    {
        if (a < b)
            return index(b, a);
        else
            return (a + columnIndex(b));
    }

    void print(){
        for (auto ii = 0u; ii < n_sites; ++ii){
            for (auto jj = ii; jj < n_sites; ++jj){
                std::cout << index(ii,jj) << " "; 
            }
            std::cout << "\n";
        }
    }
};

int main(){
    auto L = 4u;
    const std::array<unsigned int, 3> size{L,L,L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::Lattice(geometry::bccLattice(size));
    base::CouplingsMatrixDense<3> s(lat.n_sites_());
    auto nnp = base::NNProfile(lat);
    nnp.cutoff = 0.87;
    nnp.magnitude = -1.0;
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),L*L*L);
    ep.magnitude_ = 0.0;
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    //base::addDipole(s,ep);
    base::addNN(s,nnp);
    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      100;
    vmp.thermalization_ =   100000;
    vmp.n_sites_ = lat.n_sites_();
    #if 1
    std::cout << "Heisen\n";
    auto mc = montecarlo::HeatBath<3>(vmp, s);
    #endif

    #if 0
    auto v1 = geometry::Euclidean<3>({0.0, 0.0, 1.0});
    auto v2 = geometry::Euclidean<3>({0.0, 0.0, -1.0});
    auto baseset = std::vector<geometry::Euclidean<3>>({v1,v2});
    auto anivecs = montecarlo::AnisotropyVectors<3>();
    for (auto ii = 0u; ii < vmp.n_sites_; ++ii)
        anivecs.push_back(baseset);
    montecarlo::ClockModelManager<3> mc(s,anivecs);
    mc.beta_ = 1.0;
    std::cout << "Ising\n";
    #endif
    #if 0
    auto T = {  10.0, 6.0, 4.0, 3.0, 2.5,
                2.0, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9,
                0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.475, 
                0.45, 0.4, 0.35, 0.325, 0.3, 0.28, 0.26, 0.24, 0.22,
                0.20, 0.18, 0.16, 0.14, 0.12, 0.1, 0.09, 0.08,
                0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};
    #endif
    #if 1
    auto T = {0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001, 0.00000001};
    #endif
    for (auto t : T){
        auto beta = 1.0/t;
        mc.beta_ = beta;
        mc.update_(vmp.thermalization_);
        auto u = std::vector<double>(vmp.measurement_);
        auto m = std::vector<double>(vmp.measurement_);
        for (auto ii = 0u; ii < vmp.measurement_; ++ii){
            mc.update_();
            u[ii] = mc.energy_();
            m[ii] = geometry::norm(mc.spinsum_());
        }

        auto dat = montecarlo::defaultBootstrap(u, 2);
        auto energy = dat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto heatcap = dat.varianceEstimate / 
            (t * t * static_cast<double>(lat.n_sites_()));
        auto mat = montecarlo::defaultBootstrap(m, 2);
        auto magnetization = mat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi = mat.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        std::cout << t << " " << energy << " " << heatcap << " ";
        std::cout << magnetization << " " << chi << "\n";
    }
    
    std::cout << "-----------------\n";
    std::stringstream spinstring;
    mc.printSpins_(spinstring,lat);
    std::cout << spinstring.str();
}