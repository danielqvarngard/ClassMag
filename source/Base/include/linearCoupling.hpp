#ifndef CLASSMAG_BASE_LINEARCOUPLING_HPP
#define CLASSMAG_BASE_LINEARCOUPLING_HPP

#include <iostream>
#include <stdexcept>

#include "Geometry/include/matrix.hpp"
#include "spinStructure.hpp"
#include "nearestNeighbor.hpp"
#include "rkky.hpp"
#include "dipole.hpp"

namespace classmag::base{
    template<unsigned int dim>
    class LinearCouplings{
        public:
        virtual geometry::Euclidean<dim> field(
            const unsigned int site, 
            const SpinStructure<dim>& spins) const
        {
            auto result = geometry::Euclidean<dim>();
            result.fill(0.0);
            return result;
        }

        virtual inline unsigned int get_n() const {
            return 0u;
        }

        template<unsigned int latDim>
        virtual void addNN(const base::NNProfile<latDim>& nnp){

        }
        
        protected:
        LinearCouplings()
        {

        }
    };

    template<typename T, unsigned int dim>
    class CouplingsDense : public LinearCouplings<dim>{
        
        public: 
        virtual geometry::Euclidean<dim> field(
            const unsigned int site,
            const SpinStructure<dim>& spins) const override
        {
            auto result = geometry::Euclidean<dim>();
            result.fill(0.0);
            for (auto ii = 0u; ii < n_sites; ++ii){
                result += coupling(site,ii) * spins[ii];
            }
            return result;
        }

        

        inline void add(
            const unsigned int a, 
            const unsigned int b,
            const T& x)
        {
            couplingvalues[index(a,b)] += x;
        }

        inline T coupling(const unsigned int a, const unsigned int b) const
        {
            return couplingvalues[index(a,b)];
        };

        virtual inline unsigned int get_n() const override{
            return n_sites;
        }

        protected:

        std::vector<T> couplingvalues;
        CouplingsDense(){

        };

        unsigned int n_sites;

        private:
        inline unsigned int columnIndex(const unsigned int a) const
        {
            return (a * n_sites - ((a + 1u) * a)/2u);
        }
        
        inline unsigned int index(const unsigned int a, const unsigned int b) const 
        {
            if (a < b)
                return index(b, a);
            else
                return (a + columnIndex(b));
        }
    };
    
    template<unsigned int dim>
    class CouplingsMatrixDense : public CouplingsDense<geometry::Matrix<dim, dim>, dim>
    {
        public:
        CouplingsMatrixDense(){

        };

        CouplingsMatrixDense(const unsigned int n_sites)
        {
            this->n_sites = n_sites;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            this->couplingvalues.resize(n_elems);
            for (auto ii = 0u; ii < n_elems; ++ii){
                this->couplingvalues[ii] = 0.0 * geometry::eye<3>();
            }
        }

        template<unsigned int latDim>
        virtual void addNN(const NNProfile<latDim>& nnp) override
        {
            for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                    if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff)
                        this->add(ii,jj, nnp.magnitude*geometry::eye<dim>());
                }
            }
        }
    };

    template<unsigned int dim>
    class CouplingScalarDense : public CouplingsDense<double, dim>
    {
        public:
        CouplingScalarDense(const unsigned int n_sites)
        {
            this->n_sites = n_sites;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            this->couplingvalues.resize(n_elems);
            for (auto ii = 0u; ii < n_elems; ++ii){
                this->couplingvalues[ii] = 0.0;
            }
        }

        template<unsigned int dimension, unsigned int spinDim>
        virtual void addNN(const NNProfile<dimension>& nnp) override
        {
            for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                    if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff)
                        this->add(ii, jj, nnp.magnitude);
                }
            }
        };
    };
}

#endif