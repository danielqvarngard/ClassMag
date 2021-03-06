#ifndef CLASSMAG_MONTECARLO_MCPROFILE_HPP
#define CLASSMAG_MONTECARLO_MCPROFILE_HPP

#include <vector>
#include "Geometry/include/lattice.hpp"
namespace classmag::montecarlo{
    struct VectorModel_Profile{
        public:
        unsigned int thermalization_ = 100000;
        unsigned int measurement_ = 100000;
        unsigned int skips_ = 1u;
        unsigned int overrelax_ = 0u;
        unsigned int n_sites_ = 0u;
        unsigned int seed_ = 0u;
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