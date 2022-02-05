#include "include/messengers.hpp"
#include <iostream>

namespace classmag::parallelism {
    Hub::Hub()
    {
        int world_size;
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        listeners_ = world_size - 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    int Hub::gatherDoubles_(std::vector<double> &target, int tag){
        target.resize(listeners_);
        
        auto entry = 0u;
        for (int ii = 0; ii < rank_ - 1; ++ii){
            MPI_Recv(&target[entry], 1, MPI_DOUBLE, ii, tag, comm_, &status_);
            ++entry;
        }

        for (int ii = rank_ + 1; ii < listeners_ + 1; ++ii){
            MPI_Recv(&target[entry], 1, MPI_DOUBLE, ii, tag, comm_, &status_);
            ++entry;
        }

        return 0;
    }

    int Hub::gatherDoubles_(ArrayMessage &target, int tag){
        target.data_.resize(listeners_);
        
        auto entry = 0u;
        for (auto ii = 0; ii < rank_ - 1; ++ii){
            target.data_[entry].resize(target.messageLength_);
            MPI_Recv(
                &target.data_[entry][0], 
                target.messageLength_, 
                MPI_DOUBLE, 
                ii, 
                tag, 
                comm_, 
                &status_);
            ++entry;
        }

        for (auto ii = rank_ + 1; ii < listeners_ + 1; ++ii){
            target.data_[entry].resize(target.messageLength_);
            MPI_Recv(
                &target.data_[entry][0], 
                target.messageLength_, 
                MPI_DOUBLE, 
                ii, 
                tag, 
                comm_, 
                &status_);
            ++entry;
        }

        return 0;
    };

    int Hub::scatterDoubles_(const std::vector<double> &source, int tag){
//        if (source.size() != listeners_){
//            return 1;
//        }
        
        auto entry = 0u;
        for (int ii = 0; ii < rank_ - 1; ++ii){
            MPI_Send(&source[entry], 1, MPI_DOUBLE, ii, tag, comm_);
            ++entry;
        }

        for (int ii = rank_ + 1; ii < listeners_ + 1; ++ii){
            MPI_Send(&source[entry], 1, MPI_DOUBLE, ii, tag, comm_);
            ++entry;
        }

        return 0;
    }

    Listener::Listener(int hubIndex):
    hubIndex_(hubIndex)
    {

    }

    int Listener::getDouble_(double &target, int tag){
        MPI_Recv(&target, 1, MPI_DOUBLE, hubIndex_, tag, comm_, &status_);
        return 0;
    }

    int Listener::sendDouble_(const double source, int tag){
        MPI_Send(&source, 1, MPI_DOUBLE, hubIndex_, tag, comm_);
        return 0;
    }

    int Listener::sendDouble_(const std::vector<double> &source, int tag){
        for(auto d : source)
            sendDouble_(d,tag);
        return 0;
    }
}