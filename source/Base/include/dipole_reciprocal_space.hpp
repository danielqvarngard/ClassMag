#ifndef CLASSMAG_BASE_DIPOLE_RECIPROCAL_SPACE_HPP
#define CLASSMAG_BASE_DIPOLE_RECIPROCAL_SPACE_HPP

#include "Geometry/include/euclidean.hpp"
#include "Geometry/include/matrix.hpp"
#include "numerics.hpp"

namespace classmag::base
{
    geometry::Matrix<3,3> dipole_reciprocal_space_sum
    (
        const geometry::Euclidean<3> displacement,
        const double alpha,
        const unsigned int mirrors
    );


    class DipoleRecSumCalculator
    {
        public:
        DipoleRecSumCalculator(const double a, const int m):
        alpha(a), mirrors(m)
        {
            
        };

        geometry::Matrix<3,3> compute_reciprocal_space_coupling
        (
            const geometry::Euclidean<3>& displacement
        );

        private: 
        const double alpha;
        const int mirrors;
    };

    class DipoleRecSumTerm
    {
        public:
        DipoleRecSumTerm(const geometry::Euclidean<3>& d, const double a):
            displacement(d), alpha(a)
        {

        };

        geometry::Matrix<3,3> compute_term(const geometry::Euclidean<3>& n) const;

        private:
        const geometry::Euclidean<3> displacement;
        const double alpha;

        double isotropic_factor(double n_squared) const;

        inline double trig_diag(const geometry::Euclidean<3>& n) const noexcept;

        inline double trig_offdiag
        (
            const geometry::Euclidean<3>& n, 
            const unsigned int index_a, 
            const unsigned int index_b
        ) const noexcept;
    };

    

}

#endif