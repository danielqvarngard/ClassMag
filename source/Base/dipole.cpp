#include "include/dipole.hpp"

namespace classmag::base{
    geometry::Matrix<3,3> ewaldRec(
        geometry::Euclidean<3> r, 
        geometry::Euclidean<3> k, 
        double alpha){
        auto x = geometry::norm(k);
        return 4.0 * pi() * cos(k*r) * exp(-x*x/(4*alpha)) * geometry::extprod(k,k);
    }

    double ewaldRealB(double r, double alpha){
        auto result = erfc(sqrt(alpha)*r)/(r*r*r);
        result += 2.0*sqrt(alpha/pi())*exp(-alpha*r*r)/(r*r);
        return result;
    }

    double ewaldRealB(geometry::Euclidean<3> r, double alpha){
        return ewaldRealB(geometry::norm(r), alpha);
    }

    double ewaldRealC(double r, double alpha){
        auto result = 3.0 * erfc(sqrt(alpha)*r)/(pow(r,5.0));
        result += 2.0*sqrt(alpha/pi()) * (2.0 * alpha + 3/(r*r)) * exp(-alpha*r*r)/(r*r);
        return result;
    }

    double ewaldRealC(geometry::Euclidean<3> r, double alpha){
        return ewaldRealC(geometry::norm(r), alpha);
    }
}