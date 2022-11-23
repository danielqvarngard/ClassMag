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

    auto J_O = 0.0;
    auto s_O

    std::string dir = "../out/cw_";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L); 
    std::string filename = dir + label + ".txt";
    std::ofstream ofp(filename);


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
    auto c_cw = 3.0/12.0; // Curie-Weiss coefficient for a twelve-site Ising
    auto c_cw_isotropic = 1.0/(3.0* 12.0); // The same but for Heisenberg

    ofp << "Dipole: \nIsing, Heisenberg\n";

    auto sum_ising_dipole = 0.0;
    auto sum_heisen_dipole = 0.0;
    auto magnitudes = base::compute_xas0_couplings(base::RareEarthIon::Terbium, 1.0);
    for (auto ii = 0u; ii < 12u; ++ii)
    {

        for (auto jj = 0u; jj < lat.n_sites_(); ++jj)
        {
            if (jj != ii)
            {
                auto folding_index = jj % 12;
                auto dipole_matrix = magnitudes.dipole_coupling * base::dipole_fishy_matrix(ii, jj, ep);
                sum_ising_dipole += pow(anivecs[ii][0][2],2.0) * dipole_matrix[2][2] * pow(anivecs[folding_index][0][2],2.0);
                sum_heisen_dipole += dipole_matrix[2][2];
            }
        }
    }

    ofp << sum_ising_dipole << "\t" << sum_heisen_dipole << "\nRKKY:\nnu, Ising, Heisenberg\n";

    for (auto nu : nu_list)
    {
        auto sum_ising = 0.0;
        auto sum_heisen = 0.0;
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

                    sum_ising += pow(anivecs[ii][0][2],2.0) * rkky_scalar * pow(anivecs[folding_index][0][2], 2.0);
                    sum_heisen += rkky_scalar;
                }
            }
        }
        
        sum_ising *= c_cw;
        sum_heisen *= c_cw_isotropic;

        ofp << nu << "\t" << sum_ising << "\t" << sum_heisen << "\n";
    }
}