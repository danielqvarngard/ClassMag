#ifndef CLASSMAG_BASE_DIPOLE_HPP
#define CLASSMAG_BASE_DIPOLE_HPP

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "Geometry/include/lattice.hpp"
#include "Geometry/include/reciprocal.hpp"
#include "numerics.hpp"
#include "couplingLookup.hpp"
#include "linearCoupling.hpp"

namespace classmag::base{

    struct DipoleProfile{
        public:
        DipoleProfile(geometry::Lattice<3> lat):
            lattice_(lat)
            {

        };
        geometry::Lattice<3> lattice_;
        double alpha_ = 1.0;
        double magnitude_ = 1.0;
        unsigned int realMirrors_ = 0u;
        unsigned int recMirrors_ = 0u;
    };

    geometry::Euclidean<3> unscaledDipole(
        geometry::Euclidean<3> &s1, 
        geometry::Euclidean<3> &s2, 
        geometry::Euclidean<3> &r);

    double dipoleMagnitude(double moment, double latticeConstant);
    double dipoleMagnitude(double moment);

    double optimAlpha(
        const unsigned int N, 
        const unsigned int V,
        const double tauR, 
        const double tauF);

    double optimAlpha(
        const unsigned int N, 
        const unsigned int V);

    geometry::Matrix<3,3> ewaldRec(
        const geometry::Euclidean<3> &r,
        const geometry::Euclidean<3> &k, 
        const DipoleProfile &ep);
    double ewaldRealB(const double r, const DipoleProfile &ep);
    double ewaldRealB(const geometry::Euclidean<3> &r, const DipoleProfile &ep);
    geometry::Matrix<3,3> ewaldRealC(
        const geometry::Euclidean<3> &r, 
        const DipoleProfile &ep);

    geometry::Matrix<3,3> ewaldSelf(const DipoleProfile &ep);  

    geometry::Matrix<3,3> dipoleMatrix(
        const unsigned int site1, 
        const unsigned int site2,
        const DipoleProfile &ep);

    void addDipole(CouplingsMatrixDense<3> &target, const DipoleProfile& ep);
}

#endif