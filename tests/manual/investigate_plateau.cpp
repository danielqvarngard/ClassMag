#include <iostream>
#include <vector>
#include <array>
#include <string>

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/scancmd.hpp"
#include "Base/include/dipole.hpp"
#include "Base/include/linearCoupling.hpp"
#include "Base/include/default_orderparameters.hpp"
#include "Base/include/numerics.hpp"
#include "MonteCarlo/include/ClockModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"
#include "MonteCarlo/include/HeatBath.hpp"


using namespace classmag;

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

double volume(const geometry::Lattice<3> &lat){
    auto bravais = lat.getBravais_();
    auto size = lat.getSize_();
    auto V = static_cast<double>(size[0]*size[1]*size[2]);
    V *= abs(bravais[0] * geometry::cross(bravais[1],bravais[2]));
    return V;
}


int main(int argc, char* argv[])
{
    auto L = 2u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }

    auto Tl = 0.01;
    auto s_Tl = fileio::get_cmd_flag_str(argc,argv,"-Tl");
    if (!s_Tl.empty()){
        Tl = std::stoi(s_Tl);
    }

    auto Th = 10.0;
    auto s_Th = fileio::get_cmd_flag_str(argc,argv,"-Th");
    if (!s_Th.empty()){
        Th = std::stoi(s_Th);
    }

    auto n_T = 20;
    auto s_n_T = fileio::get_cmd_flag_str(argc,argv,"-n_T");
    if (!s_n_T.empty()){
        n_T = std::stoi(s_n_T);
    }

    auto T = base::linspace(Th, Tl, n_T);

    auto skip = 1;
    auto s_skip = fileio::get_cmd_flag_str(argc,argv,"-skip");
    if (!s_skip.empty()){
        skip = std::stoi(s_skip);
    }

    auto therm = 1000000;
    auto s_therm = fileio::get_cmd_flag_str(argc,argv,"-therm");
    if (!s_therm.empty()){
        therm = std::stoi(s_therm);
    }

    auto meas = 1000000;
    auto s_meas = fileio::get_cmd_flag_str(argc,argv,"-meas");
    if (!s_meas.empty()){
        meas = std::stoi(s_meas);
    }

    

    auto size = std::array<unsigned int,3>{L, L, L};
    auto lat = geometry::chas0(size);

    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      1000000;
    vmp.thermalization_ =   1000000;
    auto n_resamples = 1;
    vmp.n_sites_ = lat.n_sites_();
    
    
    
    #if 0
    auto lat = geometry::pyrochlore(size);
    auto axes = pyrochlore_ising_axes();
    auto v_set = pyrochlore_ising_axes();
    auto anivecs = montecarlo::EasyAxisVectors<3>();

    for (auto ii = 0u; ii < lat.n_sites_(); ++ii){
        auto v_p =  1.0 * v_set[ii % v_set.size()];
        auto v_m = -1.0 * v_set[ii % v_set.size()];
        auto baseset = std::vector<geometry::Euclidean<3>>({v_p, v_m});
        anivecs.push_back(baseset);
    }

    montecarlo::ClockModelManager<3> mc(s,anivecs);
    #endif
    auto ep = base::DipoleProfile(lat);
    //ep.alpha_ = base::optimAlpha(lat.n_sites_(),volume(lat));
    ep.alpha_ = 1.0;
    ep.magnitude_ = 1.0;
    auto nmax = 10u;
    ep.realMirrors_ = nmax;
    ep.recMirrors_ = nmax;
    ep.length_ = static_cast<double>(L);
    
    base::CouplingsMatrixDense s(lat.n_sites_());
    s.addDipole_spherical(ep);

    montecarlo::HeatBath<3> mc(vmp,s);

    auto delim = "-----------------\n";
    std::string dir = "../out/";
    std::string label = fileio::datestamp();
    std::string filename = dir + "plateau_"+ label + ".mcout";
    std::ofstream ofp(filename);

    auto impatience = 0.0;
    for (auto t : T)
    {
        auto beta = 1.0/t;
        mc.beta_ = beta;
        mc.update_(therm);
        auto m = std::vector<double>(meas);

        for (auto ii = 0; ii < meas; ++ii)
        {
            mc.update_(skip);
            auto v = mc.spinsum_();
            m[ii] = geometry::norm(v);
        }

        auto mat = montecarlo::defaultBootstrap(m, 10);
        auto mag = mat.meanEstimate / static_cast<double>(lat.n_sites_());
        auto chi = mat.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << t << " " << mag << " " << chi << std::endl;
        impatience = impatience + 1.0;
        std::cout << impatience/static_cast<double>(T.size()) << std::endl;
    }
    ofp.close();
    
}