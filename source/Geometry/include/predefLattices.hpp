#ifndef CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP
#define CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP

#include "lattice.hpp"

namespace classmag::geometry{
    std::vector<Euclidean<3>> icosahedralCluster();

    Lattice<2> squareLattice(const std::array<unsigned int,2> &systemSize);

    template<unsigned int dim>
    Lattice<dim> cubicLattice(
        const std::array<unsigned int, dim> &systemSize){
        std::array<Euclidean<dim>,dim> bravais;
        for (unsigned int ii = 0; ii < dim; ++ii){
            bravais[ii].fill(0.0);
            bravais[ii][ii] = 1.0;
        }

        auto lattice = Lattice<dim>(bravais,systemSize);
        return lattice;
    }
}

#endif