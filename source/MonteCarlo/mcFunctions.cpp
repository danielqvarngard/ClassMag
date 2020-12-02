#include "include/mcFunctions.hpp"

namespace classmag::interactions{
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
        double logArgument = 1.0 + z * (exponential(2.0 * fieldStrength) - 1.0);
        return (1.0/fieldStrength) * logarithm(logArgument) - 1.0;
    }
}