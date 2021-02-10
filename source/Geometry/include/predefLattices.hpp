#ifndef CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP
#define CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP

#include "lattice.hpp"

namespace classmag::geometry{
    inline double tsaiRadius();
    inline double goldenRatio();
    std::vector<Euclidean<3>> icosahedralCluster();

    SubLattice<3> bccLattice(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> has0(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> has100(const std::array<unsigned int, 3> &systemSize);

    template<unsigned int dim>
    SubLattice<dim> cubicLattice(
        const std::array<unsigned int, dim> &systemSize){
        std::array<Euclidean<dim>,dim> bravais;
        for (unsigned int ii = 0; ii < dim; ++ii){
            bravais[ii].fill(0.0);
            bravais[ii][ii] = 1.0;
        }

        auto lattice = SubLattice<dim>(bravais,systemSize);
        return lattice;
    }
}

#endif