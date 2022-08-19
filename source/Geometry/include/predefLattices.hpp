#ifndef CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP
#define CLASSMAG_GEOMETRY_PREDEFLATTICES_HPP

#include "lattice.hpp"

namespace classmag::geometry{
    inline double tsaiRadius();
    inline double goldenRatio();
    std::vector<Euclidean<3>> icosahedralCluster();
    std::vector<Euclidean<3>> icosahedralCluster(const Euclidean<3> &displacement);

    SubLattice<3> bccLattice(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> has0(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> has100(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> chas0(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> chas0_bipartite(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> chas100(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> pyrochlore(const std::array<unsigned int, 3> &systemSize);
    Lattice<3> gaa(const std::array<unsigned int, 3> &system_size);

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