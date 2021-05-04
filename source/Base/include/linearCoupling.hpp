#ifndef CLASSMAG_BASE_LINEARCOUPLING_HPP
#define CLASSMAG_BASE_LINEARCOUPLING_HPP

#include "Geometry/include/matrix.hpp"
#include "spinStructure.hpp"

namespace classmag::base{
    template<unsigned int dim>
    class LinearCouplings{
        public:
        virtual geometry::Euclidean<dim> field(
            const unsigned int site, 
            const SpinStructure<dim>& spins) const
        {
            
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
        
        inline T coupling(const unsigned int a, const unsigned int b) const
        {
            return couplingvalues[index(a,b)];
        };
        
        unsigned int n_sites;
        std::vector<T> couplingvalues;

        protected:
        CouplingsDense(){

        };
    };
    
    template<unsigned int dim>
    class CouplingsMatrixDense : public CouplingsDense<geometry::Matrix<dim, dim>, dim>
    {
        // ISK
    };

    template<unsigned int dim>
    class CouplingScalarDense : public CouplingsDense<double, dim>
    {

    };
}

#endif