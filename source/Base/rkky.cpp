#include "include/rkky.hpp"

namespace classmag::base{
    double rkky_value(const double k_F, const double r){
        auto x = 2.0*k_F*r;
        auto coupling  = -(x*cos(x) - sin(x))/(pow(x,4.0));
        return coupling;
    };

    double rkky_value(const double k_F, const double r, const double free_path){
        auto coupling  = rkky_value(k_F, r);
        coupling *= exp(-r/free_path);
        return coupling;
    };

    double xas0_k_f(const double nu, const double a)
    {
        auto N = 176.0;
        double rho = N * nu * pow(a, -3.0);
        
        double k_F = pow(3.0 * rho * pow(M_PI, 2.0),1.0/3.0);
        return k_F;
    };

    double xas0_k_f(const double nu)
    {
        return xas0_k_f(nu, 14.7e-10);
    };

    std::function<double(const unsigned int, const unsigned int)>
        rkkyInteraction(
            const double kf,
            const geometry::Lattice<3> lattice,
            const double cutoff){
        

        auto interaction = [kf, lattice, cutoff](
            const unsigned int site1, 
            const unsigned int site2){
            
            auto r2 = lattice.squareDistance_(site1,site2);
            if (lattice.squareDistance_(site1,site2)  < cutoff*cutoff && (site1 != site2)){
                auto x = 2.0*kf*sqrt(r2);
                auto coupling  = -(x*cos(x) - sin(x))/(pow(x,4.0));
                return coupling;
            }
            else{
                return 0.0;
            }
        };

        return interaction;

    }


    double rkky_scalar_mirrored
    (
        unsigned int site1, 
        unsigned int site2, 
        const RKKYProfile& rp
    ){
        
    }

}