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

    Lattice<3> gaa(const std::array<unsigned int, 3> &system_size){
        std::array<Euclidean<3>,3> bravais;
        for (auto v : bravais){
            v.fill(0.0);
        }

        bravais[0][0] = 18.78470;
        bravais[1][1] = 23.82080;
        bravais[2][2] = 5.30100;

        auto sublattice = SubLattice<3>(bravais, system_size);

        std::vector<Euclidean<3>> sites({
            {2.341701,   15.162654,   1.328378},
            {16.443000,  8.658147,    3.972623},
            {7.050649,   8.658147,    3.978878},
            {11.734051,  15.162654,   1.322122},
            {16.443000,  3.252254,    3.972623},
            {2.341701,   20.568547,   1.328378},
            {11.734051,  20.568547,   1.322122},
            {7.050649,   3.252254,    3.978878},
            {2.358231,   3.153874,    1.360290},
            {16.426469,  20.666927,   3.940710},
            {7.034119,   20.666927,   4.010790},
            {11.750582,  3.153874,    1.290210},
            {16.426469,  15.064275,   3.940710},
            {2.358231,   8.756526,    1.360290},
            {11.750582,  8.756526,    1.290210},
            {7.034119,   15.064275,   4.010790}
        });

        sublattice.decorate_(sites);
        auto lattice = Lattice(sublattice);
        return lattice;
    }

    Lattice<3> lihof4(const std::array<unsigned int, 3> &systemSize){
        auto sublattice = cubicLattice<3>(systemSize);
        std::vector<Euclidean<3>> sites({
            {0,      0,      0},
            {0,      0.5000, 0.5000},
            {0.5000, 0,      0.5000},
            {0.5000, 0.5000, 0}
        });
        sublattice.decorate_(sites);
        auto lattice = Lattice(sublattice);
        return lattice;
    }
}