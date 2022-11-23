#include <vector>
#include <iostream>
#include <cmath>

#include "Geometry/include/predefLattices.hpp"
#include "Base/include/numerics.hpp"
#include "Base/include/rkky.hpp"
#include "Base/include/dipole_spherical.hpp"
#include "Base/include/linearCoupling.hpp"
#include "Base/include/couplings_kelvin.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/scancmd.hpp"
#include "MonteCarlo/include/PredefAnisotropyVectors.hpp"
#include "MonteCarlo/include/postProcessing.hpp"

using namespace classmag;

#if 0
geometry::matrix<3,3> bare_dipole(geometry::Euclidean<3> r)
{
    auto result = -1.0*pow(geometry::norm);
}
#endif

int main(int argc, char* argv[])
{
    auto s_j = fileio::get_cmd_flag_str(argc,argv, "-n_points");
    auto n_points = 1u;
    if (!s_j.empty()){
        n_points = std::stod(s_j);
    }
    
    double nu_lower = 1.5;
    auto s_k_lower = fileio::get_cmd_flag_str(argc,argv,"-nu_l");
    if (!s_k_lower.empty()){
        nu_lower = std::stod(s_k_lower);
    }

    double nu_higher = 1.5;
    auto s_k = fileio::get_cmd_flag_str(argc,argv,"-nu_h");
    if (!s_k.empty()){
        nu_higher = std::stod(s_k);
    }

    auto L = 1u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }

    auto j_nn = 0.0;
    auto s_nn = fileio::get_cmd_flag_str(argc,argv,"-J_NN");
    if (!s_nn.empty()){
        j_nn = std::stoi(s_nn);
    }

    std::string dir = "../out/frIsing_";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L); 
    std::string filename = dir + label + ".txt";
    std::ofstream ofp(filename);

    ofp << "J_NN = " << j_nn << "\n_______________\n";

    const std::array<unsigned int, 3> size{L,L,L};

    auto lat = geometry::chas0_bipartite(size);

    auto nu_list = base::linspace(nu_lower, nu_higher, n_points);

    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),L*L*L);
    ep.magnitude_ = 1.0;
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    ep.length_ = static_cast<double>(L);

    auto anivecs = montecarlo::compute_tsai_easyaxes();
    #if 0
    {
        auto magnitudes = base::compute_xas0_couplings(base::RareEarthIon::Terbium, 1.0);
        j_nn *= magnitudes.rkky_coupling;
    }
    #endif

    double j = 6.0;
    double c_sk = pow(pow(2.0*j + 1.0, 4.0) - 1.0, 0.5)/12.0;
    c_sk *= 1.0/12.0;
    for (auto nu : nu_list)
    {
        auto sum = 0.0;
        auto k_F = base::xas0_k_f(nu, 1.0);
        auto magnitudes = base::compute_xas0_couplings(base::RareEarthIon::Terbium, nu);
        for (auto ii = 0u; ii < 12u; ++ii)
        {

            for (auto jj = 0u; jj < lat.n_sites_(); ++jj)
            {
                if (jj != ii)
                {
                    auto folding_index = jj % 12;
                    auto r = std::sqrt(lat.squareDistance_(ii,jj));
                    auto rkky_scalar = magnitudes.rkky_coupling * base::rkky_value(k_F, r);
                    rkky_scalar *= 8.0 * k_F * k_F * k_F;
                    if (r < 0.387)
                    {
                        rkky_scalar += j_nn;
                    }

                    sum += pow(anivecs[ii][0] * (rkky_scalar * anivecs[folding_index][0]), 2.0);
                }
            }
        }

        ofp << nu << "\t" << sqrt(sum) * c_sk << "\n";
    }
}