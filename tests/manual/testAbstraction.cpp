#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <mpi.h>

#include "Geometry/include/predefLattices.hpp"
#include "Base/include/nearestNeighbor.hpp"
#include "Base/include/dipole.hpp"
#include "Base/include/rkky.hpp"
#include "Base/include/linearCoupling.hpp"
#include "Base/include/numerics.hpp"
#include "MonteCarlo/include/HeatBath.hpp"
#include "MonteCarlo/include/ClockModelManager.hpp"
#include "MonteCarlo/include/postProcessing.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/scancmd.hpp"

using namespace classmag;

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
    double k_F = 1;
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

    
    
    std::string dir = "../out/";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L) + "_" + std::to_string(world_rank); 
    std::string filename = dir + label + ".mcout";
    std::ofstream ofp(filename);
    ofp << "Wavevector = " + std::to_string(k_F) + "\n";
    ofp << "J_rkky/J_dip = " + std::to_string(j) + "\n";
    ofp << "-----------------\n";

    
    const std::array<unsigned int, 3> size{L,L,2*L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::gaa(size);
    base::CouplingsMatrixDense s(lat.n_sites_());

    auto rkkyp = base::RKKYProfile(lat);
    rkkyp.k_F = k_F;
    rkkyp.cutoff = 2.0 * static_cast<double>(L);
    rkkyp.magnitude = 8.0 * j * k_F * k_F * k_F;

    auto nnp = base::NNProfile(lat);
    nnp.cutoff = 7.01;
    nnp.magnitude = -1.0;
    auto ep = base::DipoleProfile(lat);
    ep.alpha_ = base::optimAlpha(lat.n_sites_(),L*2*L*L);
    ep.magnitude_ = 1.0;
    ep.realMirrors_ = 10u;
    ep.recMirrors_ = 10u;
    s.addDipole(ep);
    //s.addNN(nnp);
    //s.addRKKY(rkkyp);
    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      100000;
    vmp.thermalization_ =   100000;
    vmp.n_sites_ = lat.n_sites_();
    #if 1
    //ofp << "Heisen\n";
    auto mc = montecarlo::HeatBath<3>(vmp, s);
    #endif

    #if 0
    auto v1 = geometry::Euclidean<3>({0.0, 0.0, 1.0});
    auto v2 = geometry::Euclidean<3>({0.0, 0.0, -1.0});
    auto baseset = std::vector<geometry::Euclidean<3>>({v1,v2});
    auto anivecs = montecarlo::EasyAxisVectors<3>();
    for (auto ii = 0u; ii < vmp.n_sites_; ++ii)
        anivecs.push_back(baseset);
    montecarlo::ClockModelManager<3> mc(s,anivecs);
    mc.beta_ = 1.0;
    ofp << "Ising\n";
    #endif
    #if 1
    auto T = {  4.0, 3.0, 2.0, 1.5,
                1.0, 0.5, 0.25, 0.1, 
                0.05, 0.025, 0.01, 0.0075, 
                0.005, 0.0025, 0.001, 0.0005};
    #endif
    std::cout << "Starting mc loop\n";
    for (auto t : T){
        auto beta = 1.0/t;
        mc.beta_ = beta;
        mc.update_(vmp.thermalization_);
        auto u = std::vector<double>(vmp.measurement_);
        auto m_a = std::vector<double>(vmp.measurement_);
        auto m_b = std::vector<double>(vmp.measurement_);
        auto m_c = std::vector<double>(vmp.measurement_);
        for (auto ii = 0u; ii < vmp.measurement_; ++ii){
            mc.update_();
            u[ii] = mc.energy_();
            auto magnetization = mc.spinsum_();
            m_a[ii] = magnetization[0];
            m_b[ii] = magnetization[1];
            m_c[ii] = magnetization[2];
        }

        auto dat = montecarlo::defaultBootstrap(u, 2);
        auto energy = dat.meanEstimate/static_cast<double>(lat.n_sites_());
        auto heatcap = dat.varianceEstimate / 
            (t * t * static_cast<double>(lat.n_sites_()));
        
        ofp << t << " " << energy << " " << heatcap << " ";

        auto mat_a = montecarlo::defaultBootstrap(m_a, 10);
        auto magnetization_a = mat_a.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi_a = mat_a.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << magnetization_a << " " << chi_a << " ";

        auto mat_b = montecarlo::defaultBootstrap(m_b, 10);
        auto magnetization_b = mat_b.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi_b = mat_b.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << magnetization_b << " " << chi_b << " ";

        auto mat_c = montecarlo::defaultBootstrap(m_c, 10);
        auto magnetization_c = mat_c.meanEstimate/static_cast<double>(lat.n_sites_());
        auto chi_c = mat_c.varianceEstimate / 
            (t * static_cast<double>(lat.n_sites_()));
        ofp << magnetization_c << " " << chi_c << "\n";

    }
    
    ofp << "-----------------\n";
    std::stringstream spinstring;
    mc.printSpins_(spinstring,lat);
    
    ofp << spinstring.str();
    MPI_Finalize();
}