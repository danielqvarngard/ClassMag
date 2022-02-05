#include <mpi.h>

#include <vector>
#include <iostream>

int main(int argc, char* argv[]){
    MPI_Init(NULL, NULL);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size >= 2)
    {
        if (world_rank == 0){
            auto msg = std::vector<double>({13,38});
            std::cout << "Sent " << msg[0] << msg[1] << "\n";
            MPI_Send(&msg[0], 2, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        }
        if (world_rank == 1){
            std::vector<double> inbox(2);
            MPI_Status status;
            MPI_Recv(&inbox[0],2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
            std::cout << "Received " << inbox[0] << inbox[1] << "\n";
        }
    }
    MPI_Finalize();
    return 0;
}