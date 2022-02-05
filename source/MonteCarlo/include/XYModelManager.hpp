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
namespace classmag::montecarlo {

    template<unsigned int spinDim>
    class XYModelManager{
        public:

        XYModelManager(
            base::LinearBase<spinDim>& cr, 
            EasyPlaneVectors<spinDim>& av):
            couplings(cr),
            anivecs(av)
        {
            spin.resize(cr.get_n());
            theta.resize(cr.get_n());
            randomize();
        }

        double theta_distr_width = 2.0 * pi();

        virtual void update_(){
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

        double energy_() const{
            auto result = 0.0;
            for (auto ii = 0u; ii < couplings.get_n(); ++ii){
                auto field = couplings.field(ii,spin);
                result -= 0.5 * (field * spin[ii]);
            }
            return result;
        }

        geometry::Euclidean<spinDim> spinsum_(){
            auto result = geometry::Euclidean<spinDim>();
            result.fill(0.0);
            for (auto ii = 0u; ii < couplings.get_n(); ++ii){
                result += spin[ii];
            }
            return result;
        }
        void setThermRuns(const unsigned int r){
            thermRuns = r;
        }
        double beta_;

        void printSpins_() const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin[ii])
                    std::cout << s << " ";
                std::cout << "\n";
            }
        }

        void printSpins_(std::stringstream& stream) const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin[ii])
                    stream << s << " ";
            }
        }

        template<unsigned int latDim>
        void printSpins_(std::stringstream& stream, geometry::Lattice<latDim>&lattice) const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin[ii])
                    stream << s << " ";
                for (auto x : lattice.position_(ii))
                    stream << x << " ";
                stream << "\n";
            }
        }

        private:
        unsigned int thermRuns = 100000;
        const base::LinearBase<spinDim>& couplings;
        const EasyPlaneVectors<spinDim>& anivecs;
        std::mt19937 rng;
        base::SpinStructure<spinDim> spin;
        std::vector<double> theta;
        std::uniform_real_distribution<double> distr = 
            std::uniform_real_distribution<double>(0,1);

        void localupdate(const unsigned int site)
        {
            auto theta_suggested = theta_distr_width * (distr(rng) - 0.5) + theta[site];
            auto suggested = compute_xy_vector(site, theta_suggested);

            auto field = couplings.field(site, spin);
            auto deltaE = (-1.0) * field * (suggested - spin[site]);
            auto r = distr(rng);
            if (boltzmannFactor(beta_,deltaE) > r){
                spin[site] = suggested;
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
            spin[site] = suggested;
        }

        void randomize(){
            for (auto ii = 0; ii < couplings.get_n(); ++ii){
                localrandom(ii);
            }
        }
    };
}

#endif