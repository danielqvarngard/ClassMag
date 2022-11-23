#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <mpi.h>

#include "Geometry/include/predefLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/dipole_spherical.hpp"
#include "Base/include/rkky.hpp"
#include "Base/include/linearCoupling.hpp"
#include "Base/include/numerics.hpp"
#include "Base/include/default_orderparameters.hpp"
#include "MonteCarlo/include/HeatBath.hpp"
#include "MonteCarlo/include/ClockModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"
#include "MonteCarlo/include/XYModelManager.hpp"
#include "MonteCarlo/include/PredefAnisotropyVectors.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/scancmd.hpp"

using namespace classmag;


double square_vecnorm(const std::vector<double>& v){
    return std::inner_product(std::begin(v), std::end(v), std::begin(v), 0.0);
}

double magnetization_norm(const std::vector<double>& v){
    double m = 0.0;
    for (auto ii = 0; ii < 3; ++ii)
        m += v[ii] * v[ii];
    return sqrt(m);
}

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    auto s_j = fileio::get_cmd_flag_str(argc,argv, "-J");
    double j = 1.0;
    if (!s_j.empty()){
        j = std::stod(s_j);
    }

    auto s_j_ex = fileio::get_cmd_flag_str(argc,argv, "-J_ex");
    double j_ex = 0.0;
    if (!s_j_ex.empty()){
        j_ex = std::stod(s_j_ex);
    }

    auto s_delta_o = fileio::get_cmd_flag_str(argc,argv, "-delta_J_o");
    double delta_o = 0.0;
    if (!s_delta_o.empty()){
        delta_o = std::stod(s_delta_o);
    }

    auto s_delta_r = fileio::get_cmd_flag_str(argc,argv, "-delta_J_r");
    double delta_r = 0.0;
    if (!s_delta_r.empty()){
        delta_r = std::stod(s_delta_r);
    }
    
    double k_F_lower = 19.8;
    auto s_k_lower = fileio::get_cmd_flag_str(argc,argv,"-k_Fl");
    if (!s_k_lower.empty()){
        k_F_lower = std::stod(s_k_lower);
    }

    double k_F_higher = 19.8;
    auto s_k_higher = fileio::get_cmd_flag_str(argc,argv,"-k_Fh");
    if (!s_k_higher.empty()){
        k_F_higher = std::stod(s_k_higher);
    }

    auto L = 2u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }
    double k_F = k_F_lower;
    if (world_size > 1){
        auto n = static_cast<double>(world_size) - 1.0;
        auto m = static_cast<double>(world_rank);
        k_F = m * (k_F_higher - k_F_lower)/n + k_F_lower;
    }

    std::string dir = "../out/heisen_rkky_";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L) + "_" + std::to_string(world_rank); 
    std::string filename = dir + label + ".mcout";
    std::ofstream ofp(filename);
    ofp << "Wavevector = " + std::to_string(k_F) + "\n";
    ofp << "J_dip/J_rkky = " + std::to_string(j) + "\n";
    ofp << "J_nn/J_rkky = " + std::to_string(j_ex) + "\n";
    ofp << "delta_o = " + std::to_string(delta_o) + "\n";
    ofp << "delta_r = " + std::to_string(delta_r) + "\n";
    ofp << "8 k_F**3 a**3 convention\n";
    ofp << "-----------------\n";

    
    const std::array<unsigned int, 3> size{L,L,L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::chas0_bipartite(size);
    base::CouplingsMatrixDense s(lat.n_sites_());

    auto rkkyp = base::RKKYProfile(lat);
    rkkyp.k_F = k_F;
    rkkyp.cutoff = 2.0 * static_cast<double>(L);
    rkkyp.magnitude = 8.0 * k_F * k_F * k_F;

    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),L*L*L);
    ep.magnitude_ = j;
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    ep.length_ = static_cast<double>(L);

    s.addRKKY(rkkyp);
    s.addDipole_fishy(ep);

    auto nnp = base::NNProfile(lat);

    double j_nn;

    nnp.cutoff = 0.3915;
    j_nn = j_ex + delta_r;
    nnp.magnitude = j_ex + delta_r;
    s.addNN(nnp);
    nnp.cutoff = 0.387;
    nnp.magnitude = -delta_r;
    s.addNN(nnp);
    nnp.cutoff = 0.37;
    nnp.magnitude = delta_o;
    s.addNN(nnp);

    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      100000;
    vmp.thermalization_ =   100000;
    auto n_resamples = 1;
    vmp.n_sites_ = lat.n_sites_();

    auto mc = montecarlo::HeatBath<3>(vmp, s);
    auto op = base::Magnetization();
    mc.addOrderParameter_(op);
    auto staggered_boy = base::BipartiteStaggeredMagnetization(lat.partitions_());
    mc.addOrderParameter_(staggered_boy);


    auto T = {  //350.0, 340.0, 330.0, 320.0, 310.0, 300.0, 290.0, 280.0, 270.0, 260.0,
                160.0, 150.0, 140.0, 130.0, 120.0, 110.0, 100.0, 95.0, 90.0,
                //80.0, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 15.0, 10.0, 5.0,
                //4.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
                85.0, 80.0, 75.0, 60.0, 45.0, 35.0, 26.0, 20.0, 16.0, 12.0, 8.5,
                5.5, 3.0, 2.0, 1.0, 0.5, 0.1, 0.025, 0.01, 0.005, 0.0025, 0.001};

    for (auto t : T){
        auto beta = 1.0/t;
        mc.beta_ = beta;
        mc.update_(vmp.thermalization_);
        auto u = std::vector<double>(vmp.measurement_);
        auto m = std::vector<double>(vmp.measurement_);
        auto stagmag = std::vector<double>(vmp.measurement_);
        for (auto ii = 0u; ii < vmp.measurement_; ++ii){
            mc.update_();
            u[ii] = mc.energy_();
            auto ops = mc.measure_();
            m[ii] = magnetization_norm(ops);
            stagmag[ii] = ops[3];
        }

        auto dat = montecarlo::defaultBootstrap(u, n_resamples);
        auto energy = dat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto heatcap = dat.varianceEstimate / 
            (t * t * static_cast<double>(lat.n_sites_()));
        auto mat = montecarlo::defaultBootstrap(m, n_resamples);
        auto magnetization = mat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi = mat.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));

        auto mat_2 = montecarlo::defaultBootstrap(stagmag, n_resamples);
        auto stag_magnetization = mat_2.meanEstimate/static_cast<double>(lat.n_sites_());
        auto stag_chi = mat_2.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << t << " " << energy << " " << heatcap << " ";
        ofp << magnetization << " " << chi << " ";
        ofp << stag_magnetization << " " << stag_chi << "\n"; 
    }
    
    ofp << "-----------------\n";
    std::stringstream spinstring;
    mc.printSpins_(spinstring,lat);
    
    ofp << spinstring.str();
    ofp.close();
    MPI_Finalize();
}