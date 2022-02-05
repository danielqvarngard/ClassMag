#ifndef CLASSMAG_PARALLELISM_MESSENGERS_HPP
#define CLASSMAG_PARALLELISM_MESSENGERS_HPP

#include <mpi.h>
#include <vector>

#include "FileIO/include/filesystem.hpp"
namespace classmag::parallelism{

    class Hub{
        public:
        Hub();
        int scatterDoubles_(const std::vector<double> &source, int tag);
        int gatherDoubles_(std::vector<double> &target, int tag);
        int gatherDoubles_(fileio::VectorTarget &target, int tag);

        private:
        int listeners_;
        int rank_ = 0;
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
#endif