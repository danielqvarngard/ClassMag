#ifndef CLASSMAG_BASE_RKKY_HPP
#define CLASSMAG_BASE_RKKY_HPP

#include <math.h>
#include <functional>
#include <utility>

#include "Geometry/include/lattice.hpp"
#include "Base/include/spinStructure.hpp"
// #include "Base/include/linearCoupling.hpp"

namespace classmag::base{
    struct RKKYProfile{
        public:
        RKKYProfile(const geometry::PointMetric& lat):
        lattice(lat)
        {

        };
        const geometry::PointMetric& lattice;
        double cutoff = 1000;
        double magnitude = 1;
        double k_F = 1;
    };

    double rkky_value(const double k_F, const double r);

    std::function<double(const unsigned int, const unsigned int)>
        rkkyInteraction(
            const double kf,
            const geometry::Lattice<3> lattice,
            const double cutoff);
    
    #if 0
    template<typename T, unsigned int spinDim>
    void addRKKY(
        LinearCouplings<T,spinDim>& target,
        const RKKYProfile& profile
    ){

    };

    template<unsigned int spinDim>
    void addRKKY(
        CouplingScalarDense<spinDim>& target,
        const RKKYProfile& profile
    ){
        for (auto ii = 0u; ii < profile.lattice.n_sites_(); ++ii){
            for (auto jj = ii + 1; jj < profile.lattice.n_sites_(); ++jj){
                auto r2 = profile.lattice.squareDistance_(ii,jj);
                if (r2 < profile.cutoff*profile.cutoff){
                    auto strength = profile.magnitude*rkky_value(profile.k_F, sqrt(r2));
                    target.add(ii, jj, strength);
                }
            }
        }
    };

    void addRKKY(
        CouplingsMatrixDense& target,
        const RKKYProfile& profile
    ){
        for (auto ii = 0u; ii < profile.lattice.n_sites_(); ++ii){
            for (auto jj = ii + 1; jj < profile.lattice.n_sites_(); ++jj){
                auto r2 = profile.lattice.squareDistance_(ii,jj);
                if (r2 < profile.cutoff*profile.cutoff){
                    auto strength = profile.magnitude*rkky_value(profile.k_F, sqrt(r2));
                    target.add(ii,jj, strength*geometry::eye<3>());
                }
            }
        }
    };
#endif
}

#endif