#ifndef CLASSMAG_BASE_LINEARCOUPLING_HPP
#define CLASSMAG_BASE_LINEARCOUPLING_HPP

#include <iostream>
#include <stdexcept>

#include "Geometry/include/matrix.hpp"
#include "spinStructure.hpp"
#include "dipole.hpp"
#include "dipole_spherical.hpp"
#include "nearestNeighbor.hpp"
#include "rkky.hpp"

namespace classmag::base{
    template<unsigned int dim>
    class LinearBase{
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

        virtual void addNN(const NNProfile& nnp) {
            std::cout << "Base addNN called\n";
        }
        virtual void addDipole(const DipoleProfile& ep) {
            std::cout << "Base addDipole called\n";
        }

        virtual void addDipole_spherical(const DipoleProfile& ep) {
            std::cout << "Base addDipole called\n";
        }

        virtual void addDipole_fishy(const DipoleProfile& ep) {
            std::cout << "Base addDipole called\n";
        }

        virtual void addRKKY(const RKKYProfile& profile) {
            std::cout << "Base addRKKY called\n";
        }
        
        protected:
        LinearBase()
        {

        }
    };

    template<typename T, unsigned int dim>
    class LinearCouplings : public LinearBase<dim>{
        public: 
        virtual inline void add(
            const unsigned int a,
            const unsigned int b,
            const T& x
        ){

        }

        protected:
        LinearCouplings()
        {

        }
    };


    // TODO: Is the spin dimension parameter necessary?
    template<typename T, unsigned int dim>
    class CouplingsDense : public LinearCouplings<T, dim>{
        
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

        

        virtual inline void add(
            const unsigned int a, 
            const unsigned int b,
            const T& x) override
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

        virtual void resize(const unsigned int n){
            n_sites = n;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            couplingvalues.resize(n_elems);
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
    
//    template<unsigned int dim>
    class CouplingsMatrixDense : public CouplingsDense<geometry::Matrix<3, 3>, 3>
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

        void resize(const unsigned int n_sites) final{
            this->n_sites = n_sites;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            this->couplingvalues.resize(n_elems);
            for (auto ii = 0u; ii < n_elems; ++ii){
                this->couplingvalues[ii] = 0.0 * geometry::eye<3>();
            }
        }

        virtual void addNN(const NNProfile& nnp) override
        {
            for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                    if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff*nnp.cutoff)
                        add(ii,jj, nnp.magnitude*geometry::eye<3>());
                }
            }
        };

        virtual void addDipole(const DipoleProfile& ep) override 
        {
            for (auto ii = 0u; ii < ep.lattice_.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < ep.lattice_.n_sites_(); ++jj){
                    auto x = dipoleMatrix(ii,jj,ep);
                    add(ii,jj,x);
                }
            }
        }

        virtual void addDipole_spherical(const DipoleProfile& ep) override 
        {
            for (auto ii = 0u; ii < ep.lattice_.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < ep.lattice_.n_sites_(); ++jj){
                    auto x = dipole_spherical_matrix(ii,jj,ep);
                    add(ii,jj,x);
                }
            }
        }

        virtual void addDipole_fishy(const DipoleProfile& ep) override 
        {
            for (auto ii = 0u; ii < ep.lattice_.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < ep.lattice_.n_sites_(); ++jj){
                    auto x = dipole_fishy_matrix(ii,jj,ep);
                    add(ii,jj,x);
                }
            }
        }

        virtual void addRKKY(const RKKYProfile& profile) override
        {
        for (auto ii = 0u; ii < profile.lattice.n_sites_(); ++ii){
            for (auto jj = ii + 1; jj < profile.lattice.n_sites_(); ++jj){
                auto r2 = profile.lattice.squareDistance_(ii,jj);
                if (r2 < profile.cutoff*profile.cutoff){
                    auto strength = profile.magnitude*rkky_value(profile.k_F, sqrt(r2));
                    add(ii,jj, strength*geometry::eye<3>());
                }
            }
        }
    };
    };

    template<unsigned int dim>
    class CouplingScalarDense : public CouplingsDense<double, dim>
    {
        public:
        CouplingScalarDense()
        {

        };

        CouplingScalarDense(const unsigned int n_sites)
        {
            this->n_sites = n_sites;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            this->couplingvalues.resize(n_elems);
            for (auto ii = 0u; ii < n_elems; ++ii){
                this->couplingvalues[ii] = 0.0;
            }
        }

        void resize(const unsigned int n_sites) final{
            this->n_sites = n_sites;
            auto n_elems = (n_sites * (n_sites + 1))/2u;
            this->couplingvalues.resize(n_elems);
            for (auto ii = 0u; ii < n_elems; ++ii){
                this->couplingvalues[ii] = 0.0;
            }
        }

        virtual void addNN(const NNProfile& nnp) override
        {
            for (auto ii = 0u; ii < nnp.lattice.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < nnp.lattice.n_sites_(); ++jj){
                    if (nnp.lattice.squareDistance_(ii,jj) < nnp.cutoff*nnp.cutoff)
                        this->add(ii,jj, nnp.magnitude);
                }
            }
        };

        virtual void addRKKY(const RKKYProfile& profile) override
        {
            for (auto ii = 0u; ii < profile.lattice.n_sites_(); ++ii){
                for (auto jj = ii + 1; jj < profile.lattice.n_sites_(); ++jj){
                    auto r2 = profile.lattice.squareDistance_(ii,jj);
                    if (r2 < profile.cutoff*profile.cutoff){
                        auto strength = profile.magnitude*rkky_value(profile.k_F, sqrt(r2));
                        this->add(ii,jj, strength);
                    }
                }
            }
        };
    };
}

#endif