#ifndef CLASSMAG_GEOMETRY_MATRIX_HPP
#define CLASSMAG_GEOMETRY_MATRIX_HPP

#include "euclidean.hpp"

namespace classmag::geometry{

    template<unsigned int m, unsigned int n>
    class Matrix : public std::array<Euclidean<n>, m>{
        
        public:
        Matrix<m,n>& operator=(const Matrix<m,n>& A){
            for (auto ii = 0u; ii < m; ++ii)
                (*this)[ii] = A[ii];
            return *this;
        }

        Matrix<m,n>& operator+=(const Matrix<m,n>& A){
            for (auto ii = 0u; ii < m; ++ii)
                (*this)[ii] += A[ii];
            return *this;
        }
    };

    template<unsigned int m, unsigned int n>
    Euclidean<m> operator*(const Matrix<m,n>& A, const Euclidean<n>& x){
        Euclidean<m> result;
        for (auto ii = 0u; ii < m; ++ii){
            result[ii] = A[ii]*x;
        }
        return result;
    }

    template<unsigned int m, unsigned int n>
    Matrix<m,n> operator*(const double &c, const Matrix<m,n>& A){
        Matrix<m,n> result;
        for (auto ii = 0u; ii < m; ++ii){
            result[ii] = c*A[ii];
        }
        return result;
    }

    template<unsigned int m, unsigned int n>
    Matrix<m,n> operator*(const Matrix<m,n>& A, const double &c){
        Matrix<m,n> result;
        for (auto ii = 0u; ii < m; ++ii){
            result[ii] = c*A[ii];
        }
        return result;
    }

    template<unsigned int m, unsigned int n>
    Matrix<m,n> operator+(const Matrix<m,n>& A, const Matrix<m,n>& B){
        Matrix<m,n> result;
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < n; ++jj)
                result[ii][jj] = A[ii][jj] + B[ii][jj];
        }
        return result;
    }

    template<unsigned int m, unsigned int n>
    Matrix<m,n> operator+(const double &c, const Matrix<m,n>& A){
        Matrix<m,n> result;
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < n; ++jj)
                result[ii][jj] = c + A[ii][jj];
        }
        return result;
    }

    template<unsigned int m, unsigned int n>
    Matrix<m,n> operator+(const Matrix<m,n>& A, const double &c){
        Matrix<m,n> result;
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < n; ++jj)
                result[ii][jj] = c + A[ii][jj];
        }
        return result;
    }

    #if 0
    template<unsigned int m, unsigned int n, unsigned int o>
    Matrix<m,o> operator*(const Matrix<m,n> &A, const Matrix<n,o> &B){
        Matrix<m,o> result;
        for (auto ii = 0u; ii < m; )
    }
    #endif

    template<unsigned int m>
    Matrix<m,m> extprod(const Euclidean<m>& e1, const Euclidean<m>& e2){
        Matrix<m,m> result;
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < m; ++ii){
                result[ii][jj] = e1[ii]*e2[ii];
            }
        }
        return result;
    }

    template<unsigned int m>
    Matrix<m,m> eye(){
        Matrix<m,m> result;
        
        for (auto ii = 0u; ii < m; ++ii){
            for (auto jj = 0u; jj < m; ++jj){
                if (ii == jj)
                    result[ii][jj] = 1.0;
                else
                    result[ii][jj] = 0.0;
            }
        }
        return result;
    }

}

#endif