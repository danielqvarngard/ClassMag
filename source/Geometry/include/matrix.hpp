#ifndef CLASSMAG_GEOMETRY_MATRIX_HPP
#define CLASSMAG_GEOMETRY_MATRIX_HPP

#include "euclidean.hpp"

namespace classmag::geometry{

    template<unsigned int m, unsigned int n>
    class matrix : public std::array<Euclidean<n>, m>{
        
        public:
        matrix<m,n>& operator=(const matrix<m,n>& A){
            for (auto ii = 0u; ii < m; ++ii)
                (*this)[m] = A[m];
            return *this;
        }

        matrix<m,n>& operator+=(const matrix<m,n>& A){
            for (auto ii = 0u; ii < m; ++ii)
                (*this)[m] += A[m];
            return *this;
        }
    };

    template<unsigned int m, unsigned int n>
    Euclidean<m> operator*(const matrix<m,n>& A, const Euclidean<n>& x){
        Euclidean<m> result;
        for (auto ii = 0u; ii < m; ++ii){
            result[ii] = A[ii]*x;
        }
        return result;
    }

    #if 0
    template<unsigned int m, unsigned int n, unsigned int o>
    matrix<m,o> operator*(const matrix<m,n> &A, const matrix<n,o> &B){
        matrix<m,o> result;
        for (auto ii = 0u; ii < m; )
    }
    #endif

    template<unsigned int m>
    matrix<m,m> extprod(Euclidean<m>& e1, Euclidean<m>& e2){
        matrix<m,m> result;
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < m; ++ii){
                result[ii][jj] = e1[ii]*e2[ii];
            }
        }
    }

}

#endif