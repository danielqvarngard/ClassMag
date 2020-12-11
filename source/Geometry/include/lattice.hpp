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

        void append_(const std::vector<Euclidean<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                decoration_.push_back(targetDecoration[ii]);
        }

        void append_(const Euclidean<dimension> &targetDecoration){
                decoration_.push_back(targetDecoration);
        }

        unsigned int n_sites_() const{
            unsigned int n = decoration_.size();
            for (auto s : systemSize_)
                n *= s;
            return n;
        }

        std::vector<unsigned int> neighborCellSites_(
            const unsigned int referenceSite) const{
            
            auto n_decorations = decoration_.size();
            auto decorationIndex = referenceSite % n_decorations;
            std::vector<unsigned int> correspondingSiteIndices(2*dimension);
            
            auto cell = referenceSite/n_decorations;
            auto cellCoordinates = std::array<unsigned int, dimension>();

            auto wrappingNumber = 1;
            for (unsigned int ii = 0; ii < dimension; ++ii){
                cellCoordinates[ii] = (cell/wrappingNumber) % systemSize_[ii];
                wrappingNumber *= systemSize_[ii];
            }

            auto atLeftEdgeOf = 
                [&cellCoordinates, this]
                (unsigned int direction){
                if (cellCoordinates[direction] % 
                        systemSize_[direction] == 0)
                    return true;
                else
                    return false;
            };
            
            auto atRightEdgeOf = 
                [&cellCoordinates, this]
                (unsigned int direction){
                if ((cellCoordinates[direction] - 1) % 
                        systemSize_[direction] == 0)
                    return true;
                else
                    return false;
            };

            auto sideIndex = 
                [referenceSite, &wrappingNumber, this]
                (unsigned int direction, bool left){
                    auto a = (left? -1 : +1);
                return (referenceSite + a * 
                        decoration_.size() * wrappingNumber);
            };

            auto wrappedSideIndex = 
                [referenceSite, &wrappingNumber, this]
                (unsigned int direction, bool left){
                    auto a = (left? 1 : -1);
                return (referenceSite + a *
                        decoration_.size() * (wrappingNumber * 
                            (systemSize_[direction] - 1)));
            };


            wrappingNumber = 1;
            for (unsigned int direction = 0; direction < dimension; ++direction){
                auto left = true;
                auto right = !left;

                if (atLeftEdgeOf(direction)){
                    correspondingSiteIndices[2*direction] = 
                        wrappedSideIndex(direction, left);
                    correspondingSiteIndices[2*direction + 1] =
                        sideIndex(direction, right); 
                }
                else if (atRightEdgeOf(direction)){
                    correspondingSiteIndices[2*direction] = 
                        sideIndex(direction, left);
                    correspondingSiteIndices[2*direction + 1] =
                        wrappedSideIndex(direction, right);
                }
                else{
                    correspondingSiteIndices[2*direction] = 
                        sideIndex(direction, left);
                    correspondingSiteIndices[2*direction + 1] =
                        sideIndex(direction, right);
                }
                wrappingNumber *= systemSize_[direction];
            }

            return correspondingSiteIndices;
        }
    };
}

#endif