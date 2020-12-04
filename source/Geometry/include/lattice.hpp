#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP

#include "euclidean.hpp"
#include <array>
#include <utility>
#include <vector>
#include <iostream>
#include <algorithm> // *min_element

namespace classmag::geometry{

    template<unsigned int dimension>
    std::vector<Euclidean<dimension>> defaultDecoration(){
        auto e = Euclidean<dimension>();
        e.fill(0.0);
        std::vector<Euclidean<dimension>> result;
        result.resize(1);
        result[0] = e;
        return result;
    }

    

    template <unsigned int dimension>
    class Lattice{
        private:
        const std::array<Euclidean<dimension>,dimension> bravais_;
        const std::array<unsigned int,dimension> systemSize_; 
        std::vector<Euclidean<dimension>> decoration_;

        public:
        Lattice<dimension>(
            const std::array<Euclidean<dimension>,dimension> &bravais,
            const std::array<unsigned int,dimension> &systemSize):
            bravais_(bravais),
            systemSize_(systemSize),
            decoration_(defaultDecoration<dimension>())
        {
        }

        Euclidean<dimension> position_(unsigned int site) const{
            auto periodicityConstant = 1;
            auto decorationNumber = site % decoration_.size();
            const int cellNumber = site / decoration_.size();
            Euclidean<dimension> translationVector = decoration_[decorationNumber];

            for (unsigned int ii = 0; ii < dimension; ++ii){
                    
                int n = (cellNumber/periodicityConstant) % systemSize_[ii];
                translationVector += n*bravais_[ii];
                periodicityConstant *= systemSize_[ii];
            }
            return translationVector;
        }

        Euclidean<dimension> mirroredPosition_(
            unsigned int site, 
            std::array<int, dimension> &periods) const{
                auto position = position_(site);
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    position += (static_cast<int>(systemSize_[ii])*periods[ii])*bravais_[ii];
                }
                return position;
        }


        virtual double squareDistance_(unsigned int site1, unsigned int site2)const {
            std::vector<double> squareDistances(3*dimension);
            for (unsigned int index = 0; index < 3*dimension; ++index){
                std::array<int, dimension> periods;
                periods.fill(0);
                int multiplier = (index % 3) - 1;
                periods[index/3] = multiplier;
                auto v1 = position_(site1);
                auto v2 = mirroredPosition_(site2,periods);
                squareDistances[index] = (v1 - v2) * (v1 - v2);
            }

            double r = *std::min_element(squareDistances.begin(), squareDistances.end());
            return r;
        }

        void decorate_(const std::vector<Euclidean<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            decoration_.resize(n_decorations);
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                decoration_[ii] = targetDecoration[ii];
        }

        unsigned int n_sites_(){
            unsigned int n = decoration_.size();
            for (auto s : systemSize_)
                n *= s;
            return n;
        }

        std::vector<unsigned int> neighborCellSites_(
            const unsigned int referenceSite){
                
            std::vector<unsigned int> correspondingSiteIndices(2*dimension);
            


        }
    };
}

#endif