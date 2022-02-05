#include <mpi.h>

#include <iostream>

#include "Base/include/linearCoupling.hpp"
#include "FileIO/include/filesystem.hpp"
#include "FileIO/include/readCouplings.hpp"
#include "FileIO/include/readMCO.hpp"
#include "FileIO/include/scancmd.hpp"
#include "FileIO/include/readDefaultSettings.hpp"
#include "Parallelism/include/messengers.hpp"
#include "Geometry/include/predefLattices.hpp"
#include "Environments/ParallelTempering/ptProcs.hpp"
#include "MonteCarlo/include/PredefAnisotropyVectors.hpp"

using namespace classmag;

int main(int argc, char* argv[]){

    MPI_Init(NULL,NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    auto infile = fileio::get_cmd_flag_str(argc,argv, "-input");
    if (infile.empty()){
        infile = "./input/default.mcinput";
    }

    auto size = fileio::readSize<3>(infile);
    auto lattice = geometry::chas0_bipartite(size);

    auto n_sites = lattice.n_sites_();

    auto temperature_intervals = fileio::read_temperatures(infile);
    auto temperatures = fileio::convert_temperature_intervals(
        world_size - 1, temperature_intervals);
    auto betas = montecarlo::invert_vector_elements(temperatures);
    montecarlo::VectorModel_Profile mco;

    fileio::readMCO(mco, infile);
    mco.n_sites_ = lattice.n_sites_();
    if (world_rank == 0){
        auto coupling_data = fileio::ReadDefaultMonteCarloOptions(infile);
        const std::string outfile = "./test.out";
        fileio::WriteDefaultFileHeader(outfile, coupling_data);
        environments::mpiPT_hub(betas,mco.measurement_);
    }
    else {
        
        auto couplings = base::CouplingsMatrixDense(n_sites);
        fileio::readLinearInteractions<geometry::Matrix<3,3>,3,3>(couplings, infile, lattice);
        auto anivecs = montecarlo::compute_tsai_easyplanes();
        auto mc_proc = montecarlo::XYModelManager<3>(couplings,anivecs);
        mc_proc.theta_distr_width *= 0.1;
        environments::mpiPT_mc<3>(mc_proc,mco);
    }
    MPI_Finalize();
}
