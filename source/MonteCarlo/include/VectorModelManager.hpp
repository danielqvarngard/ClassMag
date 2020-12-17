#ifndef CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP
#define CLASSMAG_MONTECARLO_MONTECARLOMANAGER_HPP

#include <vector>
#include <iostream>
#include <functional>
#include <random>
#include "Geometry/include/euclidean.hpp"
#include "Interactions/include/couplingLookup.hpp"
#include "Interactions/include/orderParameter.hpp"
#include "mcFunctions.hpp"

namespace classmag::montecarlo{

    template<unsigned int spinDimension>
    class VectorModelManager{
        private:
            const unsigned int n_sites_;
            const models::CouplingLookup lookup_;
            models::SpinStructure<spinDimension> spin_;
            std::mt19937 randEngine_;
            models::OrderParameters<spinDimension> orderParameters_;

            std::normal_distribution<double> normalDistribution_;
            std::uniform_real_distribution<double> uniformDistribution_;

            virtual geometry::Euclidean<spinDimension> randomUnitVector_(){
                geometry::Euclidean<spinDimension> e;
                for (unsigned int ii = 0; ii < spinDimension; ++ii){
                    e[ii] = normalDistribution_(randEngine_);
                }
                return (e/(geometry::norm(e)));
            }

            geometry::Euclidean<spinDimension> localField_(unsigned int site) const{
                geometry::Euclidean<spinDimension> field;
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

                    auto field = localField_(ii);

                    auto delta_E = (-1.0) * field * (candidateSpin - spin_[site]);
                    
                    auto r = uniformDistribution_(randEngine_);
                    auto bf = boltzmannFactor(beta_,delta_E);
                    auto check = false;
                    if (bf > r){
                        spin_[site] = candidateSpin;
                        check = true;
                    }

                }
            }

        public:
            double beta_ = 1.0;

            VectorModelManager(
                unsigned int n_sites,
                const std::function<double(unsigned int, unsigned int)> interaction, 
                int seed):
            n_sites_(n_sites),
            lookup_(models::CouplingLookup(n_sites,interaction)),
            randEngine_(std::mt19937(seed)),
            uniformDistribution_(std::uniform_real_distribution<double>()),
            normalDistribution_(std::normal_distribution<double>())
            {
                spin_.resize(n_sites);
                for (unsigned int ii = 0; ii < n_sites; ++ii)
                    spin_[ii] = randomUnitVector_();
            }

            void addOrderParameter_(
                const models::OrderParameter<spinDimension> &op){
                orderParameters_.push_back(op);
            }

            void addOrderParameter_(
                const models::OrderParameters<spinDimension> &op){
                for (auto o : op)
                    addOrderParameter_(o);
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
                for (unsigned int ii = 0; ii < n_sites_; ++ii){
                    auto field = localField_(ii);

                    spin_[ii] = 
                        2.0*(classmag::geometry::proj(spin_[ii],field)) -
                        spin_[ii];
                }

            }

            void overRelax_(const unsigned int n_times){
                for (unsigned int ii = 0; ii < n_times; ++ii)
                    overRelax_();
            }

            void update_(unsigned int n_times, unsigned int n_overRelax){
                for (unsigned int ii = 0; ii < n_times; ++ii){
                    overRelax_(n_overRelax);
                    update_();
                }
            }
            
            double energy_() const{
                auto total = 0.0;
                for (unsigned int ii = 0; ii < n_sites_; ++ii){
                    auto energy = -1.0 * localField_(ii) * spin_[ii];
                    total += energy;
                }

                return total/2.0; // Do not count interactions twice!
            }

            std::vector<double> measure_() const{
                std::vector<double> results;
                for (auto orderParameter : orderParameters_)
                    results.push_back(orderParameter.measure_(spin_));
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
                    std::cout << spin_[site][ii] << " ";
                std::cout << "\n";
            }

            void printM_() const{
                for (unsigned int ii = 0; ii < spinDimension; ++ii){
                    auto total = 0.0;
                    for (auto s : spin_)
                        total += s[ii];
                    std::cout << total << " ";
                }
                std::cout << "\n";
            }

            void printCouplings_(unsigned int site) const{
                    for (unsigned int ii = 0; ii < n_sites_; ++ii)
                        std::cout << lookup_.coupling_(site,ii) << "\n";
            }

            void printCoordinations_() const{
                    for (unsigned int ii = 0; ii < n_sites_; ++ii){
                        auto z = 0.0;
                        for (unsigned int jj = 0; jj < n_sites_; ++jj){
                            z += lookup_.coupling_(ii,jj);
                        }
                        std::cout << z << " ";
                    }
                    std::cout << "\n";
                    
                }

            void printSpins_() const{
                for (unsigned int ii = 0; ii < n_sites_; ++ii){
                    for (auto s : spin_[ii])
                        std::cout << s << " ";
                    std::cout << "\n";
                }
            }
    };

    template<>
    void VectorModelManager<3>::updateSpin_(unsigned int site){
        auto field = localField_(site);
        auto fieldNorm = geometry::norm(field);
        auto z = uniformDistribution_(randEngine_);
        
        auto cosine = heatBathCosine(beta_*fieldNorm,z);
        /*if (cosine > 1.0)
            cosine = 1.0; */ // Uncomment if returns NaN
        auto sine = squareRoot(1 - cosine*cosine);

        auto invariantAngle = 2.0*pi()*uniformDistribution_(randEngine_);
        
        auto fieldParallel = field/fieldNorm;

        auto perpendicular1 = geometry::Euclidean<3>({field[2], 0.0, -field[0]});
        perpendicular1 *= 1.0/geometry::norm(perpendicular1);
        
        auto perpendicular2 = geometry::cross(fieldParallel,perpendicular1);
        spin_[site] = 
            std::cos(invariantAngle) * sine * perpendicular1 +
            std::sin(invariantAngle) * sine * perpendicular2 + 
            cosine * fieldParallel;
    };

}

#endif