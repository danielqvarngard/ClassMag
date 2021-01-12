#ifndef CLASSMAG_BASE_NEARESTNEIGHBOR_HPP
#define CLASSMAG_BASE_NEARESTNEIGHBOR_HPP

#include <functional>
#include "Geometry/include/lattice.hpp"


namespace classmag::base{
    template<unsigned int dimension>
    std::function<double(const unsigned int, const unsigned int)> 
        nearestNeighbor(
            const double J, 
            const geometry::Lattice<dimension> &lattice, 
            const double cutoff){
        std::function<double(const unsigned int, const unsigned int)> interaction = 
            [J, lattice, cutoff](const unsigned int site1, const unsigned int site2){
                if (lattice.squareDistance_(site1,site2) < cutoff*cutoff && site1 != site2)
                    return J;
                else
                    return 0.0;
            };
        return interaction;
    };
}

#endif