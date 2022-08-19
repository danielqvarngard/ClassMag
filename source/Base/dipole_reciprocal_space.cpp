#include "include/dipole_reciprocal_space.hpp"

namespace classmag::base
{

    geometry::Matrix<3,3> dipole_reciprocal_space_sum
    (
        const geometry::Euclidean<3> displacement,
        const double alpha,
        const unsigned int mirrors
    ){
        DipoleRecSumCalculator sum(alpha, mirrors);
        return sum.compute_reciprocal_space_coupling(displacement);
    };

    geometry::Matrix<3,3> 
        DipoleRecSumCalculator::compute_reciprocal_space_coupling
    (
        const geometry::Euclidean<3>& displacement
    ){

        auto result = 4.0 * pi()/3.0 * geometry::eye<3>();

        DipoleRecSumTerm term(displacement, alpha);
        
        auto no_zero_vector = 1u;
        for (auto ii = no_zero_vector; ii < maximum_z3plus_index(mirrors); ++ii)
        {
            auto index_set = z3plus_entry(ii, mirrors);
            auto overcounting_zeroes_correction = 1.0;
            
            geometry::Euclidean<3> n;
            for (auto ii = 0u; ii < 3; ++ii)
            {
                n[ii] = static_cast<double>(index_set[ii]);
                if (index_set[ii] == 0u)
                {
                    overcounting_zeroes_correction *= 0.5;
                }
            }

            result += overcounting_zeroes_correction * 
                term.compute_term(n);
        }

        return result;
    };

    geometry::Matrix<3,3> DipoleRecSumTerm::compute_term
    (
        const geometry::Euclidean<3>& n
    )   const
    {
        geometry::Matrix<3,3> trig_matrix = trig_diag(n) * geometry::eye<3>();
        for (auto ii = 0u; ii < 3; ++ii)
        {
            trig_matrix[ii][ii] *= n[ii] * n[ii];
        }

        for (auto ii = 1u; ii < 3; ++ii)
        {
            for (auto jj = 0u; jj < ii; ++jj)
            {
                trig_matrix[ii][jj] = n[ii] * n[jj] * trig_offdiag(n, ii, jj);
                trig_matrix[jj][ii] = trig_matrix[ii][jj];
            }
        }
        
        return isotropic_factor(n * n) * trig_matrix;
    };

    double DipoleRecSumTerm::isotropic_factor(double n_squared) const
    {
        auto decay_rate = pi() * pi()/(alpha * alpha);
        auto result = 32.0 * pi() * exp(-n_squared * decay_rate)/n_squared;
        return result;
    };

    inline double DipoleRecSumTerm::trig_diag
    (
        const geometry::Euclidean<3>& n
    )   const noexcept
    {
        return  cos(2 * pi() * n[0] * displacement[0]) *
                cos(2 * pi() * n[1] * displacement[1]) * 
                cos(2 * pi() * n[2] * displacement[2]);
    };

    inline double DipoleRecSumTerm::trig_offdiag
    (
        const geometry::Euclidean<3>& n, 
        const unsigned int index_a,
        const unsigned int index_b
    )   const noexcept 
    {
        auto index_c = 3u - (index_a  + index_b); 
        return  sin(2 * pi() * n[index_a] * displacement[index_a]) *
                sin(2 * pi() * n[index_b] * displacement[index_b]) * 
                cos(2 * pi() * n[index_c] * displacement[index_c]) * (-1.0);
    };
}