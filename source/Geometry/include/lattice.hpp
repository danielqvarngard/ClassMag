#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP

#include "euclidean.hpp"
#include <array>
#include <utility>
#include <vector>

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
            const std::vector<Euclidean<dimension>> decoration_;
        public:
            Lattice<dimension>(
                const std::array<Euclidean<dimension>,dimension> &bravais,
                const std::array<unsigned int,dimension> &systemSize):
                bravais_(bravais),
                systemSize_(systemSize),
                decoration_(defaultDecoration<dimension>())
            {
            }

            Euclidean<dimension> position_(unsigned int site) const
            {
                auto periodicityConstant = 1;
                Euclidean<dimension> translationVector;
                translationVector.fill(0.0);
                
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    int n = (site/periodicityConstant) % systemSize_[ii];
                    translationVector += n*bravais_[ii];
                    periodicityConstant *= systemSize_[ii];
                }
                return translationVector;
            }

            virtual double distanceSquared_(
                const unsigned int site1, 
                const unsigned int site2) const
            {
                auto v = (position_(site1) - position_(site2));
                double result = v*v;
                return result;

            }
    };
}

#endif