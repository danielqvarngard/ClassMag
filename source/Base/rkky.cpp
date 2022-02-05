#include "include/rkky.hpp"

namespace classmag::base{
    double rkky_value(const double k_F, const double r){
        auto x = 2.0*k_F*r;
        auto coupling  = -(x*cos(x) - sin(x))/(pow(x,4.0));
        return coupling;
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

}