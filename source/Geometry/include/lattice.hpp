#ifndef CLASSMAG_GEOMETRY_LATTICE_HPP
#define CLASSMAG_GEOMETRY_LATTICE_HPP

#include "qartesian.hpp"
#include <array>
#include <cerrno>

namespace classmag::geometry{

    template <unsigned int dimension>
    class Lattice{
        private:
            const std::array<Qartesian<dimension>,dimension> bravais_(dimension);
            const std::array<unsigned int,dimension> systemSize_(dimension); 
        public:
            Lattice<dimension>(
                const std::array<Qartesian<dimension>,dimension> &bravais,
                const std::array<unsigned int,dimension> &systemSize):
                bravais_(bravais),
                systemSize_(systemSize){
            }

            Lattice<dimension>(const Lattice<dimension> &lattice){
                delete;
            }

            ~Lattice<dimension>(){
                bravais_.clear();
                systemSize_.clear();
            }

            Qartesian<dimension> double position(unsigned int site) const{
                auto periodicityConstant = 1;
                Qartesian<dimension> translationVector = {0.0, 0.0, 0.0};
                
                for (unsigned int ii = 0; ii < dimension; ++ii){
                    unsigned int n = (site1/periodicityConstant) % systemSize_[ii];
                    translationVector += n*bravais_[ii];
                    periodicityConstant *= systemSize_[ii];
                }
                return translationVector;
            }

            virtual double distanceSquared(const unsigned int site1, const unsigned int site2) const{
                auto v = (position(site1) - position(site2));
                double result = v*v;
                return result;
            }
    };
}

#endif