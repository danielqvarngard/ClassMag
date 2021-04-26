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

    double ewaldRealB(geometry::Euclidean<3> e, double alpha){
        return ewaldRealB(geometry::norm(e), alpha);
    }

    geometry::Matrix<3,3> ewaldRealC(geometry::Euclidean<3> r, double alpha){
        auto x = norm(r);
        auto c = 3.0 * erfc(sqrt(alpha)*x)/(pow(x,5.0));
        c += 2.0*sqrt(alpha/pi()) * (2.0 * alpha + 3/(r*r)) * exp(-alpha*r*r)/(r*r);
        return c * geometry::extprod(r,r);
    }

    geometry::Matrix<3,3> dipoleMatrix(
        const unsigned int site1,
        const unsigned int site2,
        const EwaldProfile &ep){
        auto result = geometry::eye<3>();
        auto range = std::vector<geometry::Euclidean<3>>();
        
        if (site1 == site2)
            range = integerSweepExclude<3>(ep.realMirrors_);
        else
            range = integerSweepFull<3>(ep.realMirrors_);
        
        for (auto n : range){
            const auto size = ep.lattice_.getSize_();
            auto mirrorcoefficients = geometry::elementwise<unsigned int,3>(size,n);
            auto mirrorvector = geometry::linearCombination<3,3>(
                mirrorcoefficients,
                ep.lattice_.getBravais_());
            auto r = ep.lattice_.position_(site1) - 
                (ep.lattice_.position_(site2) - mirrorvector);
            result += ewaldRealB(r, ep.alpha_) * geometry::eye<3>() + 
                ewaldRealC(r, ep.alpha_);
        }

        range = integerSweepExclude<3>(ep.recMirrors_);
        auto reclattice = geometry::reciprocalBasis(ep.lattice_.getBravais_());
        auto r = ep.lattice_.position_(site1) - ep.lattice_.position_(site2);
        for (auto n : range){
            auto k = geometry::linearCombination<3,3>(n,reclattice);
            result += ewaldRec(r,k,ep.alpha_);
        }
        return result;
    }

    void addDipole(MatrixLookup<3> &targetLookup, const EwaldProfile &ep){
        auto nsq = ep.lattice_.n_sites_() * ep.lattice_.n_sites_();
        for (auto ii = 0u; ii < nsq; ++ii){
            auto site1 = ii % ep.lattice_.n_sites_();
            auto site2 = (ii - site1)/ep.lattice_.n_sites_();
            targetLookup.couplingTable_[ii] += dipoleMatrix(site1, site2, ep);
        }
    }
}