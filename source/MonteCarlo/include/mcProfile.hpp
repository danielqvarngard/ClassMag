#ifndef CLASSMAG_MONTECARLO_MCPROFILE_HPP
#define CLASSMAG_MONTECARLO_MCPROFILE_HPP

#include <vector>
#include "Geometry/include/lattice.hpp"
namespace classmag::montecarlo{
    struct VectorModel_Profile{
        public:
        unsigned int thermalization_;
        unsigned int measurement_;
        unsigned int skips_;
        unsigned int overrelax_;
        unsigned int n_sites_;
        unsigned int seed_;
        std::vector<unsigned int> partitions_ = {0u};

        template<unsigned int dimension>
        void getPartitions_(const geometry::Lattice<dimension> &l){
            auto p = l.partitions_();
            partitions_.resize(p.size());
            for (auto ii = 0u; ii < p.size(); ++ii){
                partitions_[ii] = p[ii];
            }
            partitions_.push_back(l.n_sites_());
        }
    };
}

#endif 