#include "include/dipole_real_space.hpp"

namespace classmag::base
{
    geometry::Matrix<3,3> dipole_real_space_sum
    (
        const geometry::Euclidean<3> displacement,
        const double alpha,
        const unsigned int mirrors
    ){
        auto folding_vectors = integerSweepFull<3>(mirrors);
        auto result = 0.0 * geometry::eye<3>();

        for (auto n : folding_vectors)
        {
            auto folded_displacement = displacement + n;
            result += dipole_real_space_term(folded_displacement, alpha);
        }

        return result;
    };

    geometry::Matrix<3,3> dipole_real_space_term
    (
        const geometry::Euclidean<3> folded_displacement,
        const double alpha
    ){
        double folded_distance = geometry::norm(folded_displacement);
        auto result =   dipole_real_scalar(folded_distance, alpha) *
                        geometry::eye<3>();
        
        result = result + (-1.0) * dipole_real_matrix_common_factor(folded_distance, alpha) * 
            dipole_real_ext_prod(folded_displacement);
            
        return result;

    };

    double dipole_real_scalar
    (
        const double folded_distance, 
        const double alpha
    ){
        double exponential_term = 2.0 * alpha/sqrt(pi()) * 
            exp(-alpha * alpha * folded_distance * folded_distance) / 
            (folded_distance * folded_distance);
        double erfc_term = erfc(alpha * folded_distance)/
            (folded_distance * folded_distance * folded_distance);

        return exponential_term + erfc_term;
    };

    double dipole_real_matrix_common_factor
    (
        const double folded_distance,
        const double alpha
    ){
        double exponential_terms = 
            alpha /sqrt(pi()) * 
            exp(-alpha * alpha * folded_distance * folded_distance) *
            (4.0 * alpha * alpha  + 
            6.0 /(folded_distance * folded_distance)) / 
            (folded_distance * folded_distance);

        double erfc_term = 3.0 * erfc(alpha * folded_distance)/pow(folded_distance, 5.0);

        return exponential_terms + erfc_term;
    };

    geometry::Matrix<3,3> dipole_real_ext_prod
    (
        const geometry::Euclidean<3> folded_displacement
    ){
        return geometry::extprod(folded_displacement, folded_displacement);
    };
}