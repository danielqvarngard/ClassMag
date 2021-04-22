#include "include/reciprocal.hpp"

namespace classmag::geometry{
    std::array<Euclidean<3>,3> reciprocalBasis(
        const std::array<Euclidean<3>,3> &bravais){
        std::array<Euclidean<3>,3> result;

        for (auto ii = 0u; ii < 3; ++ii){
            auto a = (ii + 1) % 3;
            auto b = (ii + 2) % 3;
            auto cr = cross(bravais[a], bravais[b]);
            auto norm_factor = 2.0 * M_PI / (bravais[ii] * cr);
            result[ii] = cr * norm_factor;
        }
    }
     
}