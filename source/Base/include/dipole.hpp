#ifndef CLASSMAG_BASE_DIPOLE_HPP
#define CLASSMAG_BASE_DIPOLE_HPP

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "Geometry/include/lattice.hpp"
#include "numerics.hpp"

namespace classmag::base{
    geometry::Euclidean<3> unscaledDipole(
        geometry::Euclidean<3> &s1, 
        geometry::Euclidean<3> &s2, 
        geometry::Euclidean<3> &r);

    double dipoleMagnitude(double moment, double latticeConstant);
    double dipoleMagnitude(double moment);

    geometry::Matrix<3,3> ewaldRec(
        geometry::Euclidean<3> r,
        geometry::Euclidean<3> k, 
        double alpha);
    double ewaldRealB(double r, double alpha);
    double ewaldRealB(geometry::Euclidean<3> r, double alpha);
    double ewaldRealC(double r, double alpha);
    double ewaldRealC(geometry::Euclidean<3> r, double alpha);

    class EwaldTable{

    };
}

#endif