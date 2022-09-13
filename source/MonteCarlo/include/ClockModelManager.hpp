#ifndef CLASSMAG_MONTECARLO_CLOCKMODELMANAGER_HPP
#define CLASSMAG_MONTECARLO_CLOCKMODELMANAGER_HPP

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
    class ClockModelManager : public MonteCarloBase<spinDim>{
        public:

        ClockModelManager(
            base::LinearBase<spinDim>& cr, 
            EasyAxisVectors<spinDim>& av):
            MonteCarloBase<spinDim>(cr.get_n()),
            couplings(cr),
            anivecs(av)
        {
            representation.resize(cr.get_n());
            randomize();
        }

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
                auto field = couplings.field(ii,this->spin_);
                result -= 0.5 * (field * this->spin_[ii]);
            }
            return result;
        }

        geometry::Euclidean<spinDim> spinsum_(){
            auto result = geometry::Euclidean<spinDim>();
            result.fill(0.0);
            for (auto ii = 0u; ii < couplings.get_n(); ++ii){
                result += this->spin_[ii];
            }
            return result;
        }
        void setThermRuns(const unsigned int r){
            thermRuns = r;
        }
        double beta_;

        void printSpins_() const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin_[ii])
                    std::cout << s << " ";
                std::cout << "\n";
            }
        }

        void printSpins_(std::stringstream& stream) const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin_[ii])
                    stream << s << " ";
            }
        }

        template<unsigned int latDim>
        void printSpins_(std::stringstream& stream, geometry::Lattice<latDim>&lattice) const{
            for (unsigned int ii = 0; ii < couplings.get_n(); ++ii){
                for (auto s : this->spin_[ii])
                    stream << s << " ";
                for (auto x : lattice.position_(ii))
                    stream << x << " ";
                stream << "\n";
            }
        }

        private:
        unsigned int thermRuns = 100000;
        const base::LinearBase<spinDim>& couplings;
        const EasyAxisVectors<spinDim>& anivecs;
        std::mt19937 rng;
        std::vector<unsigned int> representation;
        std::uniform_real_distribution<double> distr = 
            std::uniform_real_distribution<double>(0,1);

        void localupdate(const unsigned int site)
        {
            auto anindex = site % anivecs.size();
            auto basisnumber = anivecs[anindex].size();
            auto suggestedIndex = 0u;
            if (basisnumber > 2){
                auto intDistr = std::uniform_int_distribution<unsigned int>(1,basisnumber);
                suggestedIndex = (representation[site] + intDistr(rng)) % basisnumber;
            }
            else
                suggestedIndex = (representation[site] + 1) % basisnumber;
            auto suggested = anivecs[anindex][suggestedIndex];

            auto field = couplings.field(site, this->spin_);
            auto deltaE = (-1.0) * field * (suggested - this->spin_[site]);
            auto r = distr(rng);
            if (boltzmannFactor(beta_,deltaE) > r){
                this->spin_[site] = suggested;
                representation[site] = suggestedIndex;
            }
        }

        void localrandom(unsigned int site){
            auto anindex = site % anivecs.size();
            auto basisnumber = anivecs[anindex].size();
            auto intDistr = std::uniform_int_distribution<unsigned int>(0,basisnumber - 1);
            auto suggestedIndex = intDistr(rng);
            representation[site] = suggestedIndex;
            auto suggested = anivecs[anindex][suggestedIndex];
            this->spin_[site] = suggested;
        }

        void randomize(){
            for (auto ii = 0; ii < couplings.get_n(); ++ii){
                localrandom(ii);
            }
        }


    };

    template<unsigned int spinDim>
    class Diagonalizer{
        public:
        Diagonalizer(){

        }

        private:
        base::SpinStructure<spinDim> spin;
        base::LinearBase<spinDim>& couplings;
        EasyAxisVectors<spinDim>& set;
        double beta;
        double partitionfcn;
        bool partitioncomputed = false;
    };
}

#endif