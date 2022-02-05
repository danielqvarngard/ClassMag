#ifndef CLASSMAG_PREDEFANISOTROPYVECTORS_HPP
#define CLASSMAG_PREDEFANISOTROPYVECTORS_HPP

#include "AnisotropyVectors.hpp"

#include <iostream>

namespace classmag::montecarlo{
    EasyPlaneVectors<3> compute_tsai_easyplanes(){
        EasyPlaneVectors<3> result;
        result.resize(12);

        auto tau = (1.0 + sqrt(5.0))/2.0;
        auto v0 = geometry::Euclidean<3>({-1.0, tau, 0.0});
        auto e = geometry::Euclidean<3>({tau, 1.0, 0.0});
        e *= 1.0/geometry::norm(e);
        v0 *= 1.0/geometry::norm(v0);

        for (unsigned int ii = 0; ii < 12; ++ii){
            result[ii][0] = v0;
            auto v1 = geometry::cross(e,v0);
            v1 *= 1.0/norm(v1);
            result[ii][1] = v1;
            circShift(e);
            if ((ii + 1) % 3 == 0)
                e[1] *= -1.0;
            
            if (ii == 5)
                e *= -1.0;
            
            circShift(v0);
            if ((ii + 1) % 3 == 0)
                v0[1] *= -1.0;
            
            if (ii == 5)
                v0 *= -1.0;
        }

        return result;
    }

}

#endif