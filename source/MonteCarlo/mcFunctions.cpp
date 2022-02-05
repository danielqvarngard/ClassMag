#include "include/mcFunctions.hpp"
#include <math.h>
namespace classmag::montecarlo{
    double squareRoot(double x){
        return std::sqrt(x);
    }

    double pi(){
        return M_PI;
    }

    double exponential(double x){
        return exp(x);
    }

    double logarithm(double x){
        return log(x);
    }

    double boltzmannFactor(double beta, double delta_E){
        return exponential(-beta*delta_E);
    };

    double heatBathCosine(double fieldStrength, double z){
        double logArgument = 
            exponential(-2.0 * fieldStrength) - 
            z * (exponential(-2.0 * fieldStrength) - 1.0);
        auto cosine = (1.0/fieldStrength) * logarithm(logArgument) + 1.0;
        return cosine;
    }

    std::vector<double> invert_vector_elements(const std::vector<double>& v){
        std::vector<double> result(v.size());
        for (auto ii = 0u; ii < v.size(); ++ii)
            result[ii] = 1.0/v[ii];
        return result;
    }
}