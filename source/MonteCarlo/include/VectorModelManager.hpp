#ifndef CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP
#define CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP

#include <vector>
#include <functional>
#include <random>
#include "Geometry/include/euclidean.hpp"
#include "Interactions/include/couplingLookup.hpp"
#include "mcFunctions.hpp"

namespace classmag::montecarlo{

    template<unsigned int spinDimension>
    class VectorModelManager{
        private:
            const unsigned int n_sites_;
            const CouplingLookup lookup_;
            std::vector<classmag::geometry::Euclidean<spinDimension>> spin_;
            std::mt19937 randEngine_;

            std::normal_distribution<double> normalDistribution_;
            std::uniform_real_distribution<double> uniformDistribution_;

            virtual Euclidean<spinDimension> randomUnitVector_() const{
                classmag::geometry::Euclidean<spinDimension> e;
                for (unsigned int ii = 0; ii < spinDimension; ++ii){
                    e[ii] = normalDistribution_(randEngine_);
                }
                return (e/e.norm());
            }

            double localField_(unsigned int site) const{
                classmag::geometry::Euclidean<spinDimension> field;
                field.fill(0.0);
                for (unsigned int ii = 0; ii < n_sites_; ++ii){
                    double coupling = lookup_.coupling_(site,ii);
                    field += coupling*spin_[ii];
                }
                return field;
            }

            virtual void updateSpin_(unsigned int site){
                for (unsigned int ii = 0; ii < n_sites_; ++ii){

                    auto candidateSpin = randomUnitVector_();

                    auto field = localField_();

                    auto delta_E = field * (candidateSpin - spin_[site]);

                    auto r = uniformDistribution_(randEngine_);
                    if (boltzmannFactor(beta_,delta_E) > r)
                        spin_[site] = candidateSpin;
                }
            }

        public:
            double beta_ = 1.0;
            VectorModelManager(
                unsigned int n_sites,
                const std::function<double(int,int)> interaction, 
                int seed):
            n_sites_(n_sites),
            lookup_(CouplingLookup(n_sites,interaction)),
            randEngine_(std::mt19937(seed)),
            uniformDistribution_(std::uniform_real_distribution<double>(0.0,1.0)),
            normalDistribution_(std::normal_distribution<double>(0.0,1.0))
            {
                spin_.resize(n_sites);
                for (unsigned int ii = 0; ii < n_sites; ++ii)
                    spin_[ii] = randomSpin_();
            }

            void update_(){
                for (unsigned int site = 0; site < n_sites_; ++site)
                    updateSpin_(site);
            }

            void update_(unsigned int n_times){
                for (unsigned int ii = 0; ii < n_times; ++ii)
                    update_();
            }

            void overRelax_(){
                for (unsigned int ii = 0; ii < n_sites; ++ii){
                    auto field = localField_(ii);

                    spin_[ii] = 
                        2.0*(classmag::geometry::proj(spin_[ii],field) -
                        spin[ii];
                }

            }
            
            double energy_() const{
                auto total = 0.0;
                for (unsigned int ii = 0; ii < n_sites_; ++ii){
                    auto energy = localField_(ii) * spin_[ii];
                    total += energy;
                }

                return total/2.0; // Do not count interactions twice!
            }

    };

}

#endif