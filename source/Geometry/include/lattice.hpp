#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP

#include "euclidean.hpp"
#include <array>
#include <cerrno>

namespace classmag::geometry{

    template <unsigned int dimension>
    class Lattice{
        private:
            const std::array<Euclidean<dimension>,dimension> bravais_;
            const std::array<unsigned int,dimension> systemSize_; 
        public:
            Lattice<dimension>(
                const std::array<Euclidean<dimension>,dimension> &bravais,
                const std::array<unsigned int,dimension> &systemSize):
                bravais_(bravais),
                systemSize_(systemSize)
            {
            }

            Lattice<dimension>(const Lattice<dimension> &lattice)
            {
                delete this;
            }

            ~Lattice<dimension>()
            {
            }

            Euclidean<dimension> position_(unsigned int site) const
            {
                auto periodicityConstant = 1;
                Euclidean<dimension> translationVector;
                translationVector.set(0.0);
                
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


    template <unsigned int dimension>
    class DecoratedLattice : public Lattice<dimension>{
        private:
            std::vector<Euclidean<dimension>> decoration_({0.0});
        public:
        void decorate_(
            const std::vector<Euclidean<dimension>> &targetDecoration){
            decoration_ = targetDecoration;
        }
    };
}

#endif