#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP
#include <array>
#include <utility>
#include <vector>
#include <iostream>
#include <algorithm> // *min_element
#include <math.h> // floor

#include "euclidean.hpp"
#include "pointMetric.hpp"

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
    class SubLattice{
        private:
        std::vector<Euclidean<dimension>> decoration_;

        public:
        SubLattice<dimension>(
            const std::array<Euclidean<dimension>,dimension> &bravais,
            const std::array<unsigned int,dimension> &systemSize):
            bravais_(bravais),
            systemSize_(systemSize),
            decoration_(defaultDecoration<dimension>())
        {
        }

        const std::array<Euclidean<dimension>,dimension> bravais_;
        const std::array<unsigned int,dimension> systemSize_; 
        

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

        unsigned int n_decorations_() const{
            return decoration_.size();
        }

        std::array<unsigned int,dimension> cellCoordinates_(unsigned int referenceSite){
            auto wrappingNumber = 1;
            auto result = std::array<unsigned int, dimension>();
            {
                const auto cell = referenceSite/n_decorations_();
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    result[ii] = (cell/wrappingNumber) % systemSize_[ii];
                    wrappingNumber *= systemSize_[ii];
                }
            }
            return result;
        }
    };

    template <unsigned int dimension>
    class Lattice : public PointMetric{
        public:

        Lattice<dimension>(){

        }

        Lattice<dimension>(
            const std::array<Euclidean<dimension>,dimension> &bravais,
            const std::array<unsigned int,dimension> &systemSize):
            bravais_(bravais),
            systemSize_(systemSize)
        {
        }

        Lattice<dimension>(SubLattice<dimension> sl):
        bravais_(sl.bravais_),
        systemSize_(sl.systemSize_),
        subLattice_({sl})
        {
        }

        std::array<Euclidean<dimension>,dimension> getBravais_() const {
            return bravais_;
        }

        std::array<unsigned int, dimension> getSize_() const {
            return systemSize_;
        }

        Euclidean<dimension> position_(unsigned int site) const{

            auto offset = 0u;
            auto sl = 0u;
            auto cont = true;
            while (cont){
                if (sl + 1 < subLattice_.size()){
                    if (site + 1 > offset + subLattice_[sl].n_sites_()){
                        offset += subLattice_[sl].n_sites_();
                        ++sl;
                    }
                    else
                        cont = false;
                }
                else
                    cont = false;
            }
            auto translationVector = subLattice_[sl].position_(site - offset);
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

        double asymDistance_(unsigned int site1, unsigned int site2) const {
            std::vector<double> squareDistances(pow(2,dimension));

            auto vb = position_(site1);
            std::array<int, dimension> periods;
            for (auto ii = 0u; ii < pow(2,dimension); ++ii){
                for(auto jj = 0u; jj < dimension; ++jj){
                    periods[jj] = static_cast<int>(floor(ii/(pow(2,jj)))) % 2;
                }
                auto va = mirroredPosition_(site2, periods);
                squareDistances[ii] = (va - vb) * (va - vb);
            }

            double r = *std::min_element(squareDistances.begin(), squareDistances.end());
            return r;
        }

        virtual double squareDistance_(unsigned int site1, unsigned int site2) const {
            std::vector<double> squareDistances(2);
            squareDistances[0] = asymDistance_(site1, site2);
            squareDistances[1] = asymDistance_(site2, site1);

            double r = *std::min_element(squareDistances.begin(), squareDistances.end());
            return r;
        }

        void decorate_(const std::vector<SubLattice<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            subLattice_.resize(n_decorations);
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                subLattice_[ii] = targetDecoration[ii];
        }

        void append_(const std::vector<SubLattice<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                subLattice_.push_back(targetDecoration[ii]);
        }

        void append_(const SubLattice<dimension> &targetDecoration){
                subLattice_.push_back(targetDecoration);
        }

        unsigned int n_sites_() const{
            unsigned int n = 0;
            for (auto s : subLattice_)
                n += s.n_sites_();
            return n;
        }

        unsigned int n_decorations_() const{
            auto result = 0u;
            for (auto sl : subLattice_)
                result += sl.n_decorations_();
            return subLattice_.size();
        }

        std::array<unsigned int,dimension> cellCoordinates(unsigned int referenceSite){
            auto wrappingNumber = 1;
            auto result = std::array<unsigned int, dimension>();
            {
                const auto cell = referenceSite/n_decorations_();
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    result[ii] = (cell/wrappingNumber) % systemSize_[ii];
                    wrappingNumber *= systemSize_[ii];
                }
            }
            return result;
        }

        // Deprecated!
        virtual std::vector<unsigned int> neighborCellSites_( 
            const unsigned int referenceSite) const{
            
            auto n_decorations = subLattice_.size();
            auto decorationIndex = referenceSite % n_decorations;
            std::vector<unsigned int> correspondingSiteIndices(2*dimension);
            auto wrappingNumber = 1;
            auto cellCoordinates = std::array<unsigned int, dimension>();
            {
                const auto cell = referenceSite/n_decorations;

                
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    cellCoordinates[ii] = (cell/wrappingNumber) % systemSize_[ii];
                    wrappingNumber *= systemSize_[ii];
                }
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
                if ((cellCoordinates[direction] + 1) % 
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
                        subLattice_.size() * wrappingNumber);
            };

            auto wrappedSideIndex = 
                [referenceSite, &wrappingNumber, this]
                (unsigned int direction, bool left){
                    auto a = (left? 1 : -1);
                return (referenceSite + a *
                        subLattice_.size() * (wrappingNumber * 
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

        std::vector<unsigned int> partitions_() const{
            std::vector<unsigned int> result (subLattice_.size() + 1, 0u);
            for (auto ii = 0u; ii < subLattice_.size(); ++ii){
                result[ii + 1] = subLattice_[ii].n_sites_() + result[ii];
            }
            return result;
        }

        void setBravais_(const std::array<Euclidean<dimension>,dimension>&bravais){
            bravais_ = bravais;
        }

        void setSize_(const std::array<unsigned int, dimension>& size){
            systemSize_ = size;
        }

        private:
        std::array<Euclidean<dimension>,dimension> bravais_;
        std::array<unsigned int,dimension> systemSize_; 
        std::vector<SubLattice<dimension>> subLattice_;
    };

    

}

#endif