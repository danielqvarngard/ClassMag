#ifndef CLASSMAG_BASE_DIPOLE_HPP
#define CLASSMAG_BASE_DIPOLE_HPP

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/reciprocal.hpp"
#include "numerics.hpp"
#include "couplingLookup.hpp"

namespace classmag::base{

    struct EwaldProfile{
        geometry::Lattice<3> lattice_;
        double alpha_;
        double magnitude_;
        double realMirrors_;
        double recMirrors_;
    };

    geometry::Euclidean<3> unscaledDipole(
        geometry::Euclidean<3> &s1, 
        geometry::Euclidean<3> &s2, 
        geometry::Euclidean<3> &r);

    double dipoleMagnitude(double moment, double latticeConstant);
    double dipoleMagnitude(double moment);

    geometry::Matrix<3,3> ewaldRec(
        const geometry::Euclidean<3> &r,
        const geometry::Euclidean<3> &k, 
        const EwaldProfile &ep);
    double ewaldRealB(const double r, const EwaldProfile &ep);
    double ewaldRealB(const geometry::Euclidean<3> &r, const EwaldProfile &ep);
    geometry::Matrix<3,3> ewaldRealC(
        const geometry::Euclidean<3> &r, 
        const EwaldProfile &ep);

    geometry::Matrix<3,3> ewaldSelf(const EwaldProfile &ep);  

    geometry::Matrix<3,3> dipoleMatrix(
        const unsigned int site1, 
        const unsigned int site2,
        const EwaldProfile &ep);

    void addDipole(MatrixLookup<3> &targetLookup, const EwaldProfile &ep);
}

#endif