#ifndef CLASSMAG_MONTECARLO_XYMODELMANAGER_HPP
#define CLASSMAG_MONTECARLO_XYMODELMANAGER_HPP

#include <random>
#include <iostream>
#include <sstream>

#include "Geometry/include/lattice.hpp"
#include "Base/include/simulationProcess.hpp"
#include "Base/include/linearCoupling.hpp"
#include "mcFunctions.hpp"
#include "AnisotropyVectors.hpp"
#include "MonteCarloBase.hpp"
namespace classmag::montecarlo {

    template<unsigned int spinDim>
    class XYModelManager : public MonteCarloBase<spinDim>{
        public:

        XYModelManager(
            base::LinearBase<spinDim>& cr, 
            EasyPlaneVectors<spinDim>& av):
            MonteCarloBase<spinDim>(cr.get_n()),
            couplings(cr),
            anivecs(av)
        {
            theta.resize(cr.get_n());
            randomize();
        }

        double theta_distr_width = 2.0 * pi();

        virtual void update_() override{
            for (unsigned int site = 0u; site < couplings.get_n(); ++site)
                localupdate(site);
        }

        virtual void update_(unsigned int ntimes){
            for (unsigned int ii= 0u; ii < ntimes; ++ii)
                update_();
        }

        virtual void thermalize_(){
            update_(thermRuns);
        }

        virtual double energy_() const override{
            auto result = 0.0;
            for (auto ii = 0u; ii < couplings.get_n(); ++ii){
                auto field = couplings.field(ii,this->spin_);
                result -= 0.5 * (field * this->spin_[ii]);
            }
            return result;
        }

        void setThermRuns(const unsigned int r){
            thermRuns = r;
        }

        private:
        unsigned int thermRuns = 100000;
        const base::LinearBase<spinDim>& couplings;
        const EasyPlaneVectors<spinDim>& anivecs;
        std::mt19937 rng;
        std::vector<double> theta;
        std::uniform_real_distribution<double> distr = 
            std::uniform_real_distribution<double>(0,1);

        void localupdate(const unsigned int site)
        {
            auto theta_suggested = theta_distr_width * (distr(rng) - 0.5) + theta[site];
            auto suggested = compute_xy_vector(site, theta_suggested);

            auto field = couplings.field(site, this->spin_);
            auto deltaE = (-1.0) * field * (suggested - this->spin_[site]);
            auto r = distr(rng);
            if (boltzmannFactor(this->beta_,deltaE) > r){
                this->spin_[site] = suggested;
                theta[site] = theta_suggested;
            }
        }

        geometry::Euclidean<spinDim> compute_xy_vector(
            const unsigned int site, 
            const double phi
        ){
            auto anindex = site % anivecs.size();
            geometry::Euclidean<spinDim> result =   anivecs[anindex][0]*cos(phi) + 
                            anivecs[anindex][1]*sin(phi);

            return result;
        }

        void localrandom(unsigned int site){
            theta[site] = 2 * pi() * distr(rng);
            auto suggested = compute_xy_vector(site, theta[site]);
            this->spin_[site] = suggested;
        }

        void randomize(){
            for (auto ii = 0; ii < couplings.get_n(); ++ii){
                localrandom(ii);
            }
        }
    };
}

#endif