#ifndef CLASSMAG_GEOMETRY_EUCLIDEAN_HPP
#define CLASSMAG_GEOMETRY_EUCLIDEAN_HPP

#include <array>
#include <utility>
#include <math.h>

namespace classmag::geometry{

    template<unsigned int dimension>
    class Euclidean : public std::array<double, dimension>{
        public:
        void fill(const double d){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] = d;
        }
        Euclidean& operator=(const Euclidean &src){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] = src[ii];
            return *this;
        }

        Euclidean& operator+=(const Euclidean &src){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] += src[ii];
            return *this;
        }

        Euclidean& operator-=(const Euclidean &src){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] -= src[ii];
            return *this;
        }

        Euclidean& operator*=(const double d){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] *= d;
            return *this;
        }
    };

    template<unsigned int dimension>
    Euclidean<dimension> operator+(
        const Euclidean<dimension> &e1,
        const Euclidean<dimension> &e2){
        Euclidean<dimension> result;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result[ii] = e1[ii] + e2[ii];
        return result;
    }

    template<unsigned int dimension>
    Euclidean<dimension> operator-(
        const Euclidean<dimension> &e1,
        const Euclidean<dimension> &e2){
        Euclidean<dimension> result;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result[ii] = e1[ii] - e2[ii];
        return result;
    }

    template<unsigned int dimension>
    double operator*(
        const Euclidean<dimension> &e1,
        const Euclidean<dimension> &e2){
        auto result = 0.0;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result += e1[ii] * e2[ii];
        return result;
    }

    template<unsigned int dimension>
    Euclidean<dimension> operator*(
        const double d,
        const Euclidean<dimension> &e){
        Euclidean<dimension> result;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result[ii] = d * e[ii];
        return result;
    }

    template<unsigned int dimension>
    Euclidean<dimension> operator*(
        const int d,
        const Euclidean<dimension> &e){
        Euclidean<dimension> result;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result[ii] = static_cast<double>(d) * e[ii];
        return result;
    }

    template<unsigned int dimension>
    Euclidean<dimension> operator*(
        const Euclidean<dimension> &e,
        const int d){
        Euclidean<dimension> result;
        for (unsigned int ii = 0; ii < dimension; ++ii)
            result[ii] = static_cast<double>(d) * e[ii];
        return result;
    }

    template<unsigned int dimension>
    Euclidean<dimension> operator/(
        const Euclidean<dimension> e,
        const double d){
        return (1.0/d) * e;
    }

    template<unsigned int dimension>
    double norm(const Euclidean<dimension> e){
        return sqrt(e*e);
    }

    template<unsigned int dimension>
    Euclidean<dimension> proj(
        const Euclidean<dimension> &e,
        const Euclidean<dimension> &onReference){
        
        Euclidean<dimension> result;
        result = ((onReference*e)/(onReference*onReference))*onReference;
        return result;
    }
}



#endif