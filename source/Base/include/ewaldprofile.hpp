#ifndef CLASSMAG_BASE_EWALDPROFILE_HPP
#define CLASSMAG_BASE_EWALDPROFILE_HPP

namespace classmag::base{
    template<unsigned int latDim>
    struct EwaldProfile{
        public:
        EwaldProfile(geometry::Lattice<latDim> lat):
            lattice_(lat)
            {

        };
        geometry::Lattice<latDim> lattice_;
        double alpha_ = 1.0;
        double magnitude_ = 1.0;
        unsigned int realMirrors_ = 0u;
        unsigned int recMirrors_ = 0u;
    };
}

#endif