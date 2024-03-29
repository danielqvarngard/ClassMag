#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <mpi.h>

#include "Geometry/include/predefLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/dipole.hpp"
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

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    auto s_lower = fileio::get_cmd_flag_str(argc,argv, "-Jl");
    double j_lower = 1.0;
    if (!s_lower.empty()){
        j_lower = std::stod(s_lower);
    }
    auto s_higher = fileio::get_cmd_flag_str(argc,argv, "-Jh");
    double j_higher = 10.0;
    if (!s_higher.empty()){
        j_higher = std::stod(s_higher);
    }
    double k_F = 19.8;
    auto s_k = fileio::get_cmd_flag_str(argc,argv,"-k_F");
    if (!s_k.empty()){
        k_F = std::stod(s_k);
    }
    auto L = 2u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }
    double j = j_lower;
    if (world_size > 1){
        auto n = static_cast<double>(world_size) - 1.0;
        auto m = static_cast<double>(world_rank);
        j = m * (j_higher - j_lower)/n + j_lower;
    }

    
    
    std::string dir = "../out/xy";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L) + "_" + std::to_string(world_rank); 
    std::string filename = dir + label + ".mcout";
    std::ofstream ofp(filename);
    ofp << "Wavevector = " + std::to_string(k_F) + "\n";
    ofp << "J_rkky/J_dip = " + std::to_string(j) + "\n";
    ofp << "-----------------\n";

    
    const std::array<unsigned int, 3> size{L,L,L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::chas0_bipartite(size);
    base::CouplingsMatrixDense s(lat.n_sites_());

    auto rkkyp = base::RKKYProfile(lat);
    rkkyp.k_F = k_F;
    rkkyp.cutoff = 2.0 * static_cast<double>(L);
    rkkyp.magnitude = 8.0 * j * k_F * k_F * k_F;

    auto nnp = base::NNProfile(lat);
    nnp.cutoff = 0.387;
    nnp.magnitude = j;
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),L*L*L);
    ep.magnitude_ = 1.0;
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    s.addDipole(ep);
    //s.addNN(nnp);
    s.addRKKY(rkkyp);
    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      10000;
    vmp.thermalization_ =   100000;
    vmp.n_sites_ = lat.n_sites_();
    #if 1
    auto anivecs = montecarlo::compute_tsai_easyplanes();
    auto mc = montecarlo::XYModelManager<3>(s, anivecs);
    auto op = base::Magnetization();
    mc.addOrderParameter_(op);
    mc.theta_distr_width = 2.0*montecarlo::pi()/20.0;
    #endif

    #if 0
    auto v1 = geometry::Euclidean<3>({0.0, 0.0, 1.0});
    auto v2 = geometry::Euclidean<3>({0.0, 0.0, -1.0});
    auto baseset = std::vector<geometry::Euclidean<3>>({v1,v2});
    auto anivecs = montecarlo::EasyAxisVectors<3>();
    for (auto ii = 0u; ii < vmp.n_sites_; ++ii)
        anivecs.push_back(baseset);
    montecarlo::ClockModelManager<3> mc(s,anivecs);

    montecarlo::
    mc.beta_ = 1.0;
    ofp << "Ising\n";
    #endif
    #if 1
    auto T = {  20.0, 16.0, 14.0, 12.0,
                9.5, 8.0, 7.0,
                5.5, 4.5, 3.5, 3.0, 2.5,
                2.0, 1.0, 0.5, 0.1, 0.025, 0.01, 0.0075, 0.005, 0.0025, 0.001};
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
            m[ii] = square_vecnorm(mc.measure_());
        }

        auto dat = montecarlo::defaultBootstrap(u, 2);
        auto energy = dat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto heatcap = dat.varianceEstimate / 
            (t * t * static_cast<double>(lat.n_sites_()));
        auto mat = montecarlo::defaultBootstrap(m, 2);
        auto magnetization = mat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi = mat.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << t << " " << energy << " " << heatcap << " ";
        ofp << magnetization << " " << chi << "\n";
    }
    
    ofp << "-----------------\n";
    std::stringstream spinstring;
    mc.printSpins_(spinstring,lat);
    
    ofp << spinstring.str();
    MPI_Finalize();
}