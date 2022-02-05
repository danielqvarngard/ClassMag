#ifndef CLASSMAG_BASE_NEARESTNEIGHBOR_HPP
#define CLASSMAG_BASE_NEARESTNEIGHBOR_HPP

#include <functional>
#include <iostream>
//#include "linearCoupling.hpp"
#include "Geometry/include/lattice.hpp"

namespace classmag::base{

    struct NNProfile{
        public:
        NNProfile(const geometry::PointMetric& latref):
        lattice(latref)
        {

        }

        const geometry::PointMetric& lattice;
        double cutoff;
        double magnitude = 1.0;
    };

    template<unsigned int dimension>
    std::function<double(const unsigned int, const unsigned int)> 
        nearestNeighbor(
            const double J, 
            const geometry::Lattice<dimension> &lattice, 
            const double cutoff){
        std::function<double(const unsigned int, const unsigned int)> interaction = 
            [J, lattice, cutoff](const unsigned int site1, const unsigned int site2){
                if (lattice.squareDistance_(site1,site2) < cutoff*cutoff && site1 != site2)
                    return J;
                else
                    return 0.0;
            };
        return interaction;
    };

#if 0 
    template<typename T, unsigned int spinDim>
    void addNN(
        LinearCouplings<T,spinDim> &target, 
        const NNProfile& nnp){
        std::cout << "Base addNN called \n";
    }

    template<unsigned int spinDim>
    void addNN(
        CouplingsMatrixDense<spinDim>& target, 
        const NNProfile& nnp)
    {
        std::cout << "Matrix addNN called\n";
        for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
            for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff*nnp.cutoff)
                    target.add(ii,jj, nnp.magnitude*geometry::eye<spinDim>());
            }
        }
    };

    template<unsigned int spinDim>
    void addNN(
        CouplingScalarDense<spinDim>& target, 
        const NNProfile& nnp)
    {
        std::cout << "Scalar addNN called\n";
        for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
            for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff*nnp.cutoff)
                    target.add(ii, jj, nnp.magnitude);
            }
        }
    };
#endif
}

#endif