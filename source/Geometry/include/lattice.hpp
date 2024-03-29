#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP
#include <array>
#include <utility>
#include <vector>
#include <iostream>
#include <algorithm> // *min_element
#include <math.h> // floor
#include <string>

#include "euclidean.hpp"
#include "matrix.hpp"
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

    template<unsigned int dimension>
    std::vector<Euclidean<dimension>> integer_sweep_positive(unsigned int nmax){
        auto n_elems = static_cast<unsigned int>(pow(nmax + 1, dimension)) - 1u;
        auto result = std::vector<Euclidean<dimension>>(n_elems);
        for (auto ii = 0u; ii < n_elems; ++ii){
            for (auto jj = 0u; jj < dimension; ++jj){
                auto numerator = ii + 1u;
                auto denominator = static_cast<unsigned int>(pow(nmax + 1,jj));
                result[ii][jj] = (numerator/denominator) % (nmax + 1);
            }
        }
        return result;
    }

    template<unsigned int dimension>
    std::vector<Euclidean<dimension>> integer_sweep_full(unsigned int nmax){
        auto offsets = integer_sweep_positive<dimension>(2u * nmax);
        auto corner = Euclidean<dimension>();
        corner.fill(-1.0 * nmax);
        auto result = std::vector<Euclidean<dimension>>(offsets.size() + 1);
        result[0] = corner;
        for (auto ii = 0u; ii < offsets.size(); ++ii){
            result[ii + 1] = corner + offsets[ii];
        }
        return result;
    }
    
    template <unsigned int dimension>
    class SubLattice{
        public:
        SubLattice<dimension>(
            const std::array<Euclidean<dimension>,dimension> &bravais,
            const std::array<unsigned int,dimension> &systemSize):
            bravais_(bravais),
            systemSize_(systemSize),
            decoration_(defaultDecoration<dimension>())
        {
        }

        SubLattice<dimension>(){};

        std::vector<Euclidean<dimension>> decoration_;
        std::array<Euclidean<dimension>,dimension> bravais_;
        std::array<unsigned int,dimension> systemSize_; 

        void set_Bravais(const std::array<Euclidean<dimension>,dimension>& bravais){
            for (auto ii = 0u; ii < dimension; ++ii){
                bravais_[ii] = bravais[ii];
            }
        }

        void set_size(const std::array<unsigned int,dimension>& size){
            for (auto ii = 0u; ii < dimension; ++ii){
                systemSize_[ii] = size[ii];
            }
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

        void decorate_(const std::vector<Euclidean<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            decoration_.resize(n_decorations);
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                decoration_[ii] = targetDecoration[ii];
        }

        void decorate_(const std::vector<std::array<double,dimension>> & targetDecoration){
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

        Euclidean<dimension> get_period_vector(unsigned int direction) const noexcept
        {
            Euclidean<dimension> result;
            result.fill(0.0);

            if (direction < dimension)
            {
                for (auto ii = 0u; ii < dimension; ++ii)
                {
                    result += static_cast<double>(systemSize_[ii]) * bravais_[ii];
                }
            }
            return result;
        }


        #if 0
        std::vector<Euclidean<dimension>> compute_folding_vectors(unsigned int nmax)
        {
            std::vector<Euclidean<dimension>> result;
            auto sweep = integer_sweep_full<dimension>(nmax);

            Matrix<dimension, dimension> period_matrix;


        }
        #endif


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

        virtual double squareDistance_(
            const unsigned int site1, 
            const unsigned int site2) const override{
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

        void decorate_(const std::vector<std::array<double,dimension>> & targetDecoration){
            if (subLattice_.size() != 1)
                subLattice_.resize(1);
            subLattice_[0].decorate_(targetDecoration);
        }



        void append_(const std::vector<SubLattice<dimension>> & targetDecoration){
            auto n_decorations = targetDecoration.size();
            for (unsigned int ii = 0; ii < n_decorations; ++ii)
                subLattice_.push_back(targetDecoration[ii]);
        }

        void append_(const SubLattice<dimension> &targetDecoration){
                subLattice_.push_back(targetDecoration);
        }

        virtual unsigned int n_sites_() const override{
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

        void setBravais_(const std::vector<std::array<double, dimension>>& bravais){
            if (bravais.size() == dimension){
                for (auto ii = 0u; ii < dimension; ++ii)
                    bravais_[ii] = bravais[ii];
            }
        }

        void setSize_(const std::array<unsigned int, dimension>& size){
            systemSize_ = size;
        }

        auto getBravais_() -> std::array<Euclidean<dimension>,dimension>
        {
            return bravais_;
        }

        void getBravais_(const std::vector<std::array<double, dimension>>& bravais){
            if (bravais.size() == dimension){
                for (auto ii = 0u; ii < dimension; ++ii)
                    bravais_[ii] = bravais[ii];
            }
        }

        auto getSize_() -> std::array<unsigned int, dimension>{
            return systemSize_;
        }

        auto getSubLattice_() -> std::vector<SubLattice<dimension>>
        {
            return subLattice_;
        }

        std::string positionString(unsigned int site, const char delim) const{
            auto pos = position_(site);
            std::string result;
            for (auto x : pos){
                result += std::to_string(x) + delim;
            }
            return result;
        }

        private:
        std::array<Euclidean<dimension>,dimension> bravais_;
        std::array<unsigned int,dimension> systemSize_; 
        std::vector<SubLattice<dimension>> subLattice_;
    };

    

}

#endif