#ifndef CLASSMAG_PARALLELISM_MESSENGERS_HPP
#define CLASSMAG_PARALLELISM_MESSENGERS_HPP

#include <mpi.h>
#include <vector>

#include "ArrayMessage.hpp"
namespace classmag::parallelism{

    class CentralNodeMessenger{
        public:
        CentralNodeMessenger();
        int scatterDoubles_(const std::vector<double> &source, int tag);
        int gatherDoubles_(std::vector<double> &target, int tag);
        int gatherDoubles_(ArrayMessage &target, int tag);

        private:
        int listeners_;
        int rank_ = 0;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };

    class EdgeNodeMessenger{
        public:
        EdgeNodeMessenger() = default;
        EdgeNodeMessenger(int hubIndex);
        int getDouble_(double &target, int tag);
        int sendDouble_(const double source, int tag);
        int sendDouble_(const std::vector<double> &source, int tag);
        
        private:
        int hubIndex_ = 0;
        MPI_Comm comm_ = MPI_COMM_WORLD;
        MPI_Status status_;
    };
}
#endif