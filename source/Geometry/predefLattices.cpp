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

    std::vector<Euclidean<3>> icosahedralCluster(const Euclidean<3> &displacement){
        auto cluster = icosahedralCluster();

        for (auto ii = 0u; ii < cluster.size(); ++ii)
            cluster[ii] += displacement;
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

    Lattice<3> chas0(const std::array<unsigned int, 3> &systemSize){
        auto sublattice = cubicLattice<3>(systemSize);
        sublattice.decorate_(icosahedralCluster());
        auto e = Euclidean<3>({0.5, 0.5, 0.5});
        sublattice.append_(icosahedralCluster(e));
        auto lattice = Lattice(sublattice);
        return lattice;
    }

    Lattice<3> chas0_bipartite(const std::array<unsigned int, 3> &systemSize){
        auto sublatticeA = cubicLattice<3>(systemSize);
        sublatticeA.decorate_(icosahedralCluster());
        auto e = Euclidean<3>({0.5, 0.5, 0.5});
        auto sublatticeB = cubicLattice<3>(systemSize);
        sublatticeB.decorate_(icosahedralCluster(e));
        auto lattice = Lattice(sublatticeA);
        lattice.append_(sublatticeB);
        return lattice;
    }

    Lattice<3> chas100(const std::array<unsigned int, 3> &systemSize){
        auto lattice = chas0(systemSize);
        auto sublattice = cubicLattice<3>(systemSize);
        auto center1 = Euclidean<3>({0.0, 0.0, 0.0});
        auto center2 = Euclidean<3>({0.5, 0.5, 0.5}); 
        std::vector<Euclidean<3>> centers({center1, center2});
        sublattice.decorate_(centers);
        lattice.append_(sublattice);
        return lattice;
    }

    Lattice<3> pyrochlore(const std::array<unsigned int, 3> &systemSize){
        auto sublattice = cubicLattice<3>(systemSize);
        std::vector<Euclidean<3>> sites({
            {0.0, 0.0, 0.0}, 
            {0.25, 0.25, 0.0}, 
            {0.25, 0.0, 0.25},
            {0.0, 0.25, 0.25}, 
            {0.5, 0.5, 0.0}, 
            {0.75, 0.75, 0.0},
            {0.75, 0.5, 0.25}, 
            {0.5, 0.75, 0.25}, 
            {0.0, 0.5, 0.5}, 
            {0.25, 0.75, 0.5},
            {0.25, 0.5, 0.75}, 
            {0.0, 0.75, 0.75}, 
            {0.5, 0.0, 0.5}, 
            {0.75, 0.25, 0.5},
            {0.75, 0.0, 0.75}, 
            {0.5, 0.25, 0.75}});
        sublattice.decorate_(sites);
        auto lattice = Lattice(sublattice);
        return lattice;
    }
}