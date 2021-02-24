#ifndef CLASSMAG_GEOMETRY_DISORDEREDLATTICES_HPP
#define CLASSMAG_GEOMETRY_DISORDEREDLATTICES_HPP

#include <random>

#include "lattice.hpp"
#include "predefLattices.hpp"

namespace classmag::geometry{
    template<unsigned int dimension>
    SubLattice<dimension> disorderedCubic(
        const std::array<unsigned int, dimension> &systemSize,
        const double occupancy,
        const Euclidean<dimension> &displacement
    ){
        std::array<unsigned int, dimension> tempSize;
        tempSize.fill(1u);
        std::array<Euclidean<dimension>,dimension> tempBravais;
        for (auto ii = 0u; ii < dimension; ++ii){
            tempBravais[ii].fill(0.0);
            tempBravais[ii][ii] = static_cast<double>(systemSize[ii]);
        }

        auto result = SubLattice<dimension>(tempBravais, tempSize);
        auto proposeSites = cubicLattice<dimension>(systemSize);
        proposeSites.decorate_({displacement});
        std::vector<Euclidean<dimension>> acceptedSites;
        std::uniform_real_distribution<double> uniformDistribution = 
            std::uniform_real_distribution<double>(0,1);
        std::mt19937 rng_(std::random_device{}());

        for (auto ii = 0u; ii < proposeSites.n_sites_(); ++ii){
            if (uniformDistribution(rng_) < occupancy)
                acceptedSites.push_back(proposeSites.position_(ii));
        }

        result.decorate_(acceptedSites);
        return result;
    }

    template<unsigned int dimension>
    SubLattice<dimension> disorderedCubic(
        const std::array<unsigned int, dimension> &systemSize,
        const double occupancy
    ){
        Euclidean<dimension> disp;
        disp.fill(0.0);
        auto dc = disorderedCubic<dimension>(systemSize, occupancy, disp);
        return dc;
    }
}

#endif