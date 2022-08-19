#ifndef CLASSMAG_MONTECARLO_HEATBATH_HPP
#define CLASSMAG_MONTECARLO_HEATBATH_HPP

#include <iostream>
#include <sstream>
#include <random>

#include "Geometry/include/euclidean.hpp"
#include "Base/include/magnetizationData.hpp"
#include "Base/include/couplingLookup.hpp"
#include "Base/include/simulationProcess.hpp"
#include "Base/include/orderParameter.hpp"
#include "Base/include/linearCoupling.hpp"
#include "MonteCarloBase.hpp"
#include "mcProfile.hpp"
#include "mcFunctions.hpp"
namespace classmag::montecarlo{

    template<unsigned int spinDimension>
    class HeatBath : public MonteCarloBase<spinDimension>{
        public:
        HeatBath():
        base::SimulationBase<spinDimension>(),
        rng_(std::mt19937(mcp_.seed_)),
        normalDistribution_(std::normal_distribution<double>()),
        mcp_(VectorModel_Profile())
        {

        };

        HeatBath(
            base::LinearBase<spinDimension> &c):
        base::SimulationBase<spinDimension>(c.get_n()),
        mcp_(VectorModel_Profile()),
        lookup_(c),
        rng_(std::mt19937(mcp_.seed_)),
        normalDistribution_(std::normal_distribution<double>())
        {
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii)
                this->spin_[ii] = randomUnitVector_();
        }

        HeatBath(
            VectorModel_Profile &mcp,
            base::LinearBase<spinDimension> &c):
        MonteCarloBase<spinDimension>(mcp.n_sites_),
        mcp_(mcp),
        lookup_(c),
        rng_(std::mt19937(mcp.seed_)),
        normalDistribution_(std::normal_distribution<double>())
        {
            for (unsigned int ii = 0; ii < mcp.n_sites_; ++ii)
                this->spin_[ii] = randomUnitVector_();
        }

        HeatBath(const HeatBath &vm):
        base::SimulationBase<spinDimension>(vm.mcp_.n_sites_),
        mcp_(vm.mcp_),
        lookup_(vm.lookup_),
        rng_(std::mt19937(vm.mcp_.seed_)),
        normalDistribution_(std::normal_distribution<double>())
        {
            for (unsigned int ii = 0; ii < mcp_.n_sites_; ++ii)
                this->spin_[ii] = randomUnitVector_();
        };

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

        void update_(const unsigned int n_times){
            update_(n_times, mcp_.overrelax_);
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
            
        virtual double energy_() const override{
            auto total = 0.0;
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto energy = -1.0 * localField_(ii) * this->spin_[ii];
                total += energy;
            }

            return total/2.0; // Do not count interactions twice!
        }

        geometry::Euclidean<spinDimension> spinsum_() const{
            auto result = geometry::Euclidean<spinDimension>();
            result.fill(0.0);
            for (auto ii = 0u; ii < this->n_sites_; ++ii){
                result += this->spin_[ii];
            }
            return result;
        }

        std::vector<geometry::Euclidean<spinDimension>> magnetization_(){
            auto n_sublattices = mcp_.partitions_.size() - 1;
            geometry::Euclidean<spinDimension> initialvalue;
            initialvalue.fill(0.0);
            std::vector<geometry::Euclidean<spinDimension>>  result(n_sublattices, initialvalue);

            for (auto ii = 0u; ii < n_sublattices; ++ii){
                for (
                    auto site = mcp_.partitions_[ii]; 
                    site < mcp_.partitions_[ii + 1]; 
                    ++site){
                    
                    result[ii] += this->spin_[site];
                }
            }

            return result;

        }
        
        private:
        VectorModel_Profile mcp_;
        base::LinearBase<spinDimension>& lookup_;
        std::mt19937 rng_;

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
            return lookup_.field(site, this->spin_);
        }

        virtual void updateSpin_(unsigned int site){
            for (unsigned int ii = 0; ii < this->n_sites_; ++ii){
                auto candidateSpin = randomUnitVector_();
                auto field = localField_(ii);
                auto delta_E = (-1.0) * field * (candidateSpin - this->spin_[site]);
                auto r = uniformDistribution_(rng_);
                auto bf = boltzmannFactor(this->beta_,delta_E);
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
    void HeatBath<3>::updateSpin_(unsigned int site){
        auto field = localField_(site);
        auto fieldNorm = geometry::norm(field);
        auto z = uniformDistribution_(rng_);
        
        auto cosine = heatBathCosine(this->beta_*fieldNorm,z);
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