#include "include/dipole.hpp"

namespace classmag::base{
    geometry::Matrix<3,3> ewaldRec(
        const geometry::Euclidean<3> &r, 
        const geometry::Euclidean<3> &k, 
        const EwaldProfile &ep){
        
        auto size = ep.lattice_.getSize_();
        auto bravais = ep.lattice_.getBravais_();

        auto V = 1.0;
        for (unsigned int ii = 0u; ii < 3; ++ii){
            V *= size[ii];
        }
        V *= abs(bravais[0] * geometry::cross(bravais[1], bravais[2]));

        auto x = geometry::norm(k);
        return (4.0 * pi() / V) *
            (cos(k*r) / (x * x)) * 
            exp(-x*x/(4*ep.alpha_)) * geometry::extprod(k,k);
    }

    double ewaldRealB(const double r, const EwaldProfile &ep){
        auto result = erfc(sqrt(ep.alpha_)*r)/(r*r*r);
        result += 2.0*sqrt(ep.alpha_/pi())*exp(-ep.alpha_*r*r)/(r*r);
        return result;
    }

    double ewaldRealB(const geometry::Euclidean<3> &e, const EwaldProfile &ep){
        return ewaldRealB(geometry::norm(e), ep);
    }

    geometry::Matrix<3,3> ewaldRealC(const geometry::Euclidean<3> &r, const EwaldProfile &ep){
        auto x = norm(r);
        auto c = 3.0 * erfc(sqrt(ep.alpha_)*x)/(pow(x,5.0));
        c += 2.0*sqrt(ep.alpha_/pi()) * (2.0 * ep.alpha_ + 3/(r*r)) * exp(-ep.alpha_*r*r)/(r*r);
        return c * geometry::extprod(r,r);
    }

    geometry::Matrix<3,3> ewaldSelf(const EwaldProfile &ep){
        return 2.0 * pi()/3.0 * pow(ep.alpha_/pi(),1.5) * geometry::eye<3>();
    }

    geometry::Matrix<3,3> dipoleMatrix(
        const unsigned int site1,
        const unsigned int site2,
        const EwaldProfile &ep){
        auto result = 0.0*geometry::eye<3>();
        auto range = std::vector<geometry::Euclidean<3>>();

        if (site1 == site2){
            range = integerSweepExclude<3>(ep.realMirrors_);
            result += ewaldSelf(ep);
        }
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
            result += ewaldRealB(r, ep) * geometry::eye<3>() + 
                ewaldRealC(r, ep);
        }

        range = integerSweepExclude<3>(ep.recMirrors_);
        auto reclattice = geometry::reciprocalBasis(ep.lattice_.getBravais_());
        auto r = ep.lattice_.position_(site1) - ep.lattice_.position_(site2);
        for (auto n : range){
            auto k = geometry::linearCombination<3,3>(n,reclattice);
            result += ewaldRec(r,k,ep);
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