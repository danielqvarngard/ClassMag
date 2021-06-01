#include "include/dipole.hpp"
#include "include/linearCoupling.hpp"

namespace classmag::base{
    double optimAlpha(
        const unsigned int N, 
        const unsigned int V,
        const double tauR, 
        const double tauF)
    {
        auto arg = (tauR * pi() * N)/(tauF * V * V);
        return pow(arg, 1.0/6.0);
    }

    double optimAlpha(
        const unsigned int N, 
        const unsigned int V)
    {
        return optimAlpha(N, V, 1.0, 1.0);
    }

    geometry::Matrix<3,3> ewaldRec(
        const geometry::Euclidean<3> &r, 
        const geometry::Euclidean<3> &k, 
        const DipoleProfile &ep){
        
        auto size = ep.lattice_.getSize_();
        auto bravais = ep.lattice_.getBravais_();

        auto V = 1.0;
        for (unsigned int ii = 0u; ii < 3; ++ii){
            V *= size[ii];
        }
        V *= abs(bravais[0] * geometry::cross(bravais[1], bravais[2]));

        return (4.0 * pi() / V) *
            (cos(k*r) / (k * k)) * 
            exp(-(k * k)/(4 * ep.alpha_)) * geometry::extprod(k,k);
    }

    double ewaldRealB(const double r, const DipoleProfile &ep){
        auto result = erfc(sqrt(ep.alpha_)*r)/(r*r*r);
        result += 2.0*sqrt(ep.alpha_/pi())*exp(-ep.alpha_*r*r)/(r*r);
        return result;
    }

    double ewaldRealB(const geometry::Euclidean<3> &e, const DipoleProfile &ep){
        return ewaldRealB(geometry::norm(e), ep);
    }

    geometry::Matrix<3,3> ewaldRealC(const geometry::Euclidean<3> &r, const DipoleProfile &ep){
        auto x = norm(r);
        auto c = 3.0 * erfc(sqrt(ep.alpha_)*x)/(x*x*x*x*x);
        c += 2.0*sqrt(ep.alpha_/pi()) * (2.0 * ep.alpha_ + 3/(x*x)) * exp(-ep.alpha_*x*x)/(x*x);
        return c * geometry::extprod(r,r);
    }

    geometry::Matrix<3,3> ewaldSelf(const DipoleProfile &ep){
        return -2.0 * pi()/3.0 * pow(ep.alpha_/pi(),1.5) * geometry::eye<3>();
    }

    geometry::Matrix<3,3> dipoleMatrix(
        const unsigned int site1,
        const unsigned int site2,
        const DipoleProfile &ep){
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
                ((-1.0) * ewaldRealC(r, ep));
        }
        range = integerSweepExclude<3>(ep.recMirrors_);
        auto reclattice = geometry::reciprocalBasis(ep.lattice_.getBravais_());
        for (auto ii = 0u; ii < 3; ++ii){
            reclattice[ii] *= (M_PI/3.0)/ep.lattice_.getSize_()[ii]; // wtf is this prefactor
        }
        auto r = ep.lattice_.position_(site1) - ep.lattice_.position_(site2);
        for (auto n : range){
            auto k = geometry::linearCombination<3,3>(n,reclattice);
            result += ewaldRec(r,k,ep);
        }
        return (-1.0) * result;
    }

    void addDipole(CouplingsMatrixDense<3> &target, const DipoleProfile& ep){
        for (auto ii = 0u; ii < ep.lattice_.n_sites_(); ++ii){
            for (auto jj = ii; jj < ep.lattice_.n_sites_(); ++jj){
                auto x = dipoleMatrix(ii,jj,ep);
                target.add(ii,jj,x);
            }
        }
    }
}