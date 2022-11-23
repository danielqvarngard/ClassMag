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

    auto s_o_l = fileio::get_cmd_flag_str(argc,argv, "-J_o_l");
    double j_o_l = 0.0;
    if (!s_o_l.empty()){
        j_o_l = std::stod(s_o_l);
    }

    auto s_o_h = fileio::get_cmd_flag_str(argc,argv, "-J_o_h");
    double j_o_h = 0.0;
    if (!s_o_h.empty()){
        j_o_h = std::stod(s_o_h);
    }

    auto s_r_l = fileio::get_cmd_flag_str(argc,argv, "-J_r_l");
    double j_r_l = 0.0;
    if (!s_r_l.empty()){
        j_r_l = std::stod(s_r_l);
    }

    auto s_r_h = fileio::get_cmd_flag_str(argc,argv, "-J_r_h");
    double j_r_h = 0.0;
    if (!s_r_h.empty()){
        j_r_h = std::stod(s_r_h);
    }

    auto L = 2u;
    auto s_L = fileio::get_cmd_flag_str(argc,argv,"-L");
    if (!s_L.empty()){
        L = std::stoi(s_L);
    }

    

    std::string dir = "../out/heisen_nn_sweep_";
    std::string label = fileio::datestamp();
    label += "_" + std::to_string(L) + "_" + std::to_string(world_rank); 
    std::string filename = dir + label + ".mcout";
    std::ofstream ofp(filename);

    // Check if world_size is a square:

    auto grid_edge = static_cast<int>(floor(sqrt(world_size)));

    if (grid_edge * grid_edge != world_size)
    {
        ofp << "Warning: non square number of processors used\n";
    }



    double j_o = j_o_l;
    double j_r = j_r_l;
    if (world_size > 1){
        auto n_procs = world_size - 1;
        auto m = world_rank;
        auto n = static_cast<double>(grid_edge) - 1.0;
        
        auto x = static_cast<double>(m % grid_edge);
        auto y = static_cast<double>(m / grid_edge);

        j_o = x * (j_o_h - j_o_l)/n + j_o_l;
        j_r = y * (j_r_h - j_r_l)/n + j_r_l;
    }

    ofp << "j_o = " + std::to_string(j_o) + "\n";
    ofp << "j_i = " + std::to_string(j) + "\n";
    ofp << "j_r = " + std::to_string(j_r) + "\n";

    ofp << "-----------------\n";

    
    const std::array<unsigned int, 3> size{L,L,L};
    //auto lat = geometry::Lattice(geometry::cubicLattice<3>(size));
    auto lat = geometry::chas0_bipartite(size);
    base::CouplingsMatrixDense s(lat.n_sites_());

    //s.addRKKY(rkkyp);
    //s.addDipole_fishy(ep);

    auto nnp = base::NNProfile(lat);

    nnp.cutoff = 0.3915;
    auto delta_r = j_r - j;
    nnp.magnitude = j + delta_r;
    s.addNN(nnp);
    nnp.cutoff = 0.387;
    nnp.magnitude = -delta_r;
    s.addNN(nnp);
    auto delta_o = j_o - j;
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