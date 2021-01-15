#include <mpi.h>
#include <vector>

namespace classmag::parallelism{

    struct VectorTarget{
        public:
        VectorTarget(const unsigned int n):
        messageLength_(n)
        {

        }

        std::vector<std::vector<double>> data_;
        const unsigned int messageLength_;
    };
    class Hub{
        public:
        Hub();
        int scatterDoubles_(const std::vector<double> &source, int tag);
        int gatherDoubles_(std::vector<double> &target, int tag);
        int gatherDoubles_(VectorTarget &target, int tag);

        private:
        int listeners_;
        int rank_;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };

    class Listener{
        public:
        Listener() = default;
        Listener(int hubIndex);
        int getDouble_(double &target, int tag);
        int sendDouble_(const double source, int tag);
        int sendDouble_(const std::vector<double> &source, int tag);
        
        private:
        int hubIndex_ = 0;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };
}