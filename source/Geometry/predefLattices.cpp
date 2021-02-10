#include "include/predefLattices.hpp"
#include <iostream>

namespace classmag::geometry{
    inline double tsaiRadius(){
        return 0.3578;
    }

    inline double goldenRatio(){
        return (1.0 + sqrt(5))/2.0;
    }

    std::vector<Euclidean<3>> icosahedralCluster(){
        std::vector<Euclidean<3>> cluster(12);

        auto e = Euclidean<3>({goldenRatio(),1.0,0});
        e *= tsaiRadius()/norm(e);

        for (unsigned int ii = 0; ii < 12; ++ii){
            cluster[ii] = e;
            circShift(e);
            if ((ii + 1) % 3 == 0)
                e[1] *= -1.0;
            
            if (ii == 5)
                e *= -1.0;
        }
        return cluster;
    }

    SubLattice<3> bccLattice(const std::array<unsigned int, 3> &systemSize){
        std::array<Euclidean<3>,3> bravais;
        auto e0 = Euclidean<3>{0.5, 0.5, 0.5};
        for(unsigned int ii = 0; ii < 3; ++ii){
            bravais[ii] = e0;
            bravais[ii][ii] *= -1.0;
        }

        auto lattice = SubLattice<3>(bravais,systemSize);
        return lattice;
    }

    Lattice<3> has0(const std::array<unsigned int, 3> &systemSize){
        auto sublattice = bccLattice(systemSize);
        sublattice.decorate_(icosahedralCluster());
        auto lattice = Lattice(sublattice);
        return lattice;
    }

    Lattice<3> has100(const std::array<unsigned int, 3> &systemSize){
        auto lattice = has0(systemSize);
        lattice.append_(bccLattice(systemSize));
        return lattice;
    }
}