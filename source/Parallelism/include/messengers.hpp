#include <mpi.h>
#include <vector>

namespace classmag::parallelism{
    class Hub{
        public:
        Hub();
        int scatterDoubles_(const std::vector<double> &source, int tag);
        int gatherDoubles_(std::vector<double> &target, int tag);

        private:
        int listeners_;
        int rank_;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };

    class Listener{
        public:
        Listener(int hubIndex);
        int getDouble(double &target, int tag);
        int sendDouble(const double source, int tag);
        
        private:
        int hubIndex_ = 0;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };
}