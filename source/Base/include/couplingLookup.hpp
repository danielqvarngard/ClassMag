#ifndef CLASSMAG_BASE_COUPLINGLOOKUP_HPP
#define CLASSMAG_BASE_COUPLINGLOOKUP_HPP
#include <vector>
#include <functional>
#include <iostream>


#include "Geometry/include/matrix.hpp"
namespace classmag::base{
    class CouplingLookup{
        private:
            std::vector<double> couplingTable_;
            const unsigned int n_sites_;
        public:
            CouplingLookup(
                const unsigned int n_sites, 
                const std::function<double(int, int)> interaction);
            CouplingLookup(): n_sites_(0){};
            CouplingLookup(const CouplingLookup &src);
            double coupling_(unsigned int site1, unsigned int site2) const;
    };

    template<unsigned int spinDim>
    class MatrixLookup{
        public:
        MatrixLookup(const unsigned int n_sites):
        n_sites_(n_sites),
        couplingTable_(
            std::vector<geometry::matrix<spinDim,spinDim>>(
                n_sites*n_sites, 0.0*geometry::eye<spinDim>())
            )
        {};

        geometry::matrix<spinDim,spinDim> coupling_(unsigned int ii, unsigned int jj){
            return couplingTable_[n_sites_*site1 + site2];
        };

        private:
        const unsigned int n_sites_ = 0u;
        std::vector<geometry::matrix<spinDim,spinDim>> couplingTable_;
    };
}

#endif