#ifndef CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP
#define CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP

#include <iostream>
#include <random>

#include "Geometry/include/euclidean.hpp"
#include "Base/include/couplingLookup.hpp"
#include "Base/include/simulationProcess.hpp"
#include "Base/include/orderParameter.hpp"
#include "mcProfile.hpp"
#include "mcFunctions.hpp"
namespace classmag::montecarlo{

    template<unsigned int spinDimension>
    class VectorModelManager : public base::SimulationBase<spinDimension>{
        public:
        VectorModelManager(
            MC_Profile &mcp,
            const std::function<double(unsigned int, unsigned int)> interaction):
        base::SimulationBase<spinDimension>(mcp.n_sites_),
        mcp_(mcp),
        lookup_(base::CouplingLookup(mcp.n_sites_,interaction)),
        rng_(std::mt19937(mcp.seed_)),
        normalDistribution_(std::normal_distribution<double>())
        {
            for (unsigned int ii = 0; ii < mcp.n_sites_; ++ii)
                this->spin_[ii] = randomUnitVector_();
        }

        VectorModelManager(const VectorModelManager &vm) = delete;

        void addOrderParameter_(
            const base::OrderParameter<spinDimension> &op){
            orderParameters_.push_back(op);
        }

        void addOrderParameter_(
            const base::OrderParameters<spinDimension> &op){
            for (auto o : op)
                addOrderParameter_(o);
        }

        void overRelax_(){
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto field = localField_(ii);

                this->spin_[ii] = 
                    2.0*(geometry::proj(this->spin_[ii],field)) -
                    this->spin_[ii];
            }

        }

        void overRelax_(const unsigned int n_times){
            for (unsigned int ii = 0; ii < n_times; ++ii)
                overRelax_();
        }

        void update_(const unsigned int n_times, const unsigned int n_overRelax){
            for (unsigned int ii = 0; ii < n_times; ++ii){
                overRelax_(n_overRelax);
                updateLattice_();
            }
        }

        virtual void update_() override{
            update_(mcp_.skips_,mcp_.overrelax_);
        }

        virtual void thermalize_(){
            this->update_(mcp_.thermalization_, mcp_.overrelax_);
        }

        unsigned int measurements_(){
            return mcp_.measurement_;
        }
            
        double energy_() const{
            auto total = 0.0;
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto energy = -1.0 * localField_(ii) * this->spin_[ii];
                total += energy;
            }

            return total/2.0; // Do not count interactions twice!
        }

        virtual std::vector<double> measure_() const override{
            std::vector<double> results;
            for (auto orderParameter : orderParameters_)
                results.push_back(orderParameter.measure_(this->spin_));
            return results;
        }

        std::vector<std::string> measurementNames_() const{
            std::vector<std::string> results;
            for (auto measure : orderParameters_)
                results.push_back(measure.name_);
            return results;
        }

        void printSpin_(unsigned int site) const{
            for (unsigned int ii = 0; ii < spinDimension; ++ii)
                std::cout << this->spin_[site][ii] << " ";
            std::cout << "\n";
        }

        void printM_() const{
            for (unsigned int ii = 0; ii < spinDimension; ++ii){
                auto total = 0.0;
                for (auto s : this->spin_)
                    total += s[ii];
                std::cout << total << " ";
            }
            std::cout << "\n";
        }

        void printCouplings_(unsigned int site) const{
                for (unsigned int ii = 0; ii < this->n_sites_; ++ii)
                    std::cout << lookup_.coupling_(site,ii) << "\n";
        }

        void printCoordinations_() const{
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto z = 0.0;
                for (unsigned int jj = 0; jj < this->n_sites_; ++jj){
                    z += lookup_.coupling_(ii,jj);
                }
                std::cout << z << " ";
            }
            std::cout << "\n";        
        }

        void printSpins_() const{
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                for (auto s : this->spin_[ii])
                    std::cout << s << " ";
                std::cout << "\n";
            }
        }

        double beta_ = 1.0;
        
        private:
        const MC_Profile mcp_;
        base::CouplingLookup lookup_;
        std::mt19937 rng_;
        base::OrderParameters<spinDimension> orderParameters_;

        std::normal_distribution<double> normalDistribution_;
        std::uniform_real_distribution<double> uniformDistribution_ = 
            std::uniform_real_distribution<double>(0,1);

        virtual geometry::Euclidean<spinDimension> randomUnitVector_(){
            geometry::Euclidean<spinDimension> e;
            for (unsigned int ii = 0; ii < spinDimension; ++ii){
                e[ii] = normalDistribution_(rng_);
            }
            return (e/(geometry::norm(e)));
        }

        virtual geometry::Euclidean<spinDimension> localField_(unsigned int site) const{
            geometry::Euclidean<spinDimension> field;
            field.fill(0.0);
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                double coupling = lookup_.coupling_(site,ii);
                field += coupling*this->spin_[ii];
            }
            return field;
        }

        virtual void updateSpin_(unsigned int site){
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto candidateSpin = randomUnitVector_();
                auto field = localField_(ii);
                auto delta_E = (-1.0) * field * (candidateSpin - this->spin_[site]);
                auto r = uniformDistribution_(rng_);
                auto bf = boltzmannFactor(beta_,delta_E);
                auto check = false;
                if (bf > r){
                    this->spin_[site] = candidateSpin;
                    check = true;
                }

            }
        }

        void updateLattice_(){
            for (unsigned int site = 0; site < this->n_sites_; ++site)
                updateSpin_(site);
        }
    };

    template<>
    void VectorModelManager<3>::updateSpin_(unsigned int site){
        auto field = localField_(site);
        auto fieldNorm = geometry::norm(field);
        auto z = uniformDistribution_(rng_);
        
        auto cosine = heatBathCosine(beta_*fieldNorm,z);
        /*if (cosine > 1.0)
            cosine = 1.0; */ // Uncomment if returns NaN
        auto sine = squareRoot(1 - cosine*cosine);

        auto invariantAngle = 2.0*pi()*uniformDistribution_(rng_);
        
        auto fieldParallel = field/fieldNorm;

        auto perpendicular1 = geometry::Euclidean<3>({field[2], 0.0, -field[0]});
        perpendicular1 *= 1.0/geometry::norm(perpendicular1);
        
        auto perpendicular2 = geometry::cross(fieldParallel,perpendicular1);
        this->spin_[site] = 
            std::cos(invariantAngle) * sine * perpendicular1 +
            std::sin(invariantAngle) * sine * perpendicular2 + 
            cosine * fieldParallel;
    };

}

#endif