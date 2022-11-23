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

    auto s_j = fileio::get_cmd_flag_str(argc,argv, "-J_I");
    double j = 1.0;
    if (!s_j.empty()){
        j = std::stod(s_j);
    }
    
    double j_o_lower = -1.0;
    auto s_o_lower = fileio::get_cmd_flag_str(argc,argv,"-J_O_l");
    if (!s_o_lower.empty()){
        j_o_lower = std::stod(s_o_lower);
    }

    double j_o_higher = 1.0;
    auto s_o_higher = fileio::get_cmd_flag_str(argc,argv,"-J_O_h");
    if (!s_o_higher.empty()){
        j_o_higher = std::stod(s_o_higher);
    }

    auto L = 2u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }

    double j_o = j_o_lower;
    if (world_size > 1){
        auto n = static_cast<double>(world_size) - 1.0;
        auto m = static_cast<double>(world_rank);
        j_o = m * (j_o_higher - j_o_lower)/n + j_o_lower;
    }

    std::string dir = "../out/ising_nn";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L) + "_" + std::to_string(world_rank); 
    std::string filename = dir + label + ".mcout";
    std::ofstream ofp(filename);
    ofp << "J_I = " + std::to_string(j) + "\n";
    ofp << "J_O = " + std::to_string(j_o) + "\n";
    ofp << "-----------------\n";

    
    const std::array<unsigned int, 3> size{L,L,L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::chas0_bipartite(size);
    base::CouplingScalarDense<3> s(lat.n_sites_());

    auto nnp = base::NNProfile(lat);
    nnp.cutoff = 0.387;
    nnp.magnitude = j;
    s.addNN(nnp);
    nnp.cutoff = 0.37;
    auto delta = j_o - j;
    nnp.magnitude = delta;
    s.addNN(nnp);
    auto vmp = montecarlo::VectorModel_Profile();
    vmp.measurement_ =      50000;
    vmp.thermalization_ =   100000;
    auto n_resamples = 1;
    vmp.n_sites_ = lat.n_sites_();
    #if 1
    auto anivecs = montecarlo::compute_tsai_easyaxes();
    auto mc = montecarlo::ClockModelManager<3>(s, anivecs);
    auto op = base::Magnetization();
    mc.addOrderParameter_(op);
    auto staggered_boy = base::BipartiteStaggeredMagnetization(lat.partitions_());
    mc.addOrderParameter_(staggered_boy);
    #endif

    #if 0
    auto v1 = geometry::Euclidean<3>({0.0, 0.0, 1.0});
    auto v2 = geometry::Euclidean<3>({0.0, 0.0, -1.0});
    auto baseset = std::vector<geometry::Euclidean<3>>({v1,v2});
    auto anivecs = montecarlo::EasyAxisVectors<3>();
    for (auto ii = 0u; ii < vmp.n_sites_; ++ii)
        anivecs.push_back(baseset);
    montecarlo::ClockModelManager<3> mc(s,anivecs);

    ofp << "Ising\n";
    #endif
    #if 1
    auto T = {  20.0, 15.0, 10.0,
                5.0, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 
                0.9, 0.8, 0.7,
                0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15,
                0.14, 0.13, 0.12, 0.11,
                0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 
                0.025, 0.01, 0.005, 0.0025, 0.001};
                //4.0, 3.0, 2.0, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
                //85.0, 80.0, 75.0, 60.0, 45.0, 35.0, 26.0, 20.0, 16.0, 12.0, 8.5,
    #endif
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
            stagmag[ii] = abs(ops[3]);
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