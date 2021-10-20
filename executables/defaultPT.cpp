#include <mpi.h>

#include "Environments/ParallelTempering/ptProcs.hpp"
#include "FileIO/include/optionsReadin.hpp"

using namespace classmag;

int main(int argc, char* argv[]){
    MPI_Init(NULL,NULL);

    int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size == 1)
        return -1;
    if (world_rank == 0)
        environments::mpiPT_hub({1.0, 2.0}, 0);
    else{
        
    }
}