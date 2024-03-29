#ifndef CLASSMAG_GEOMETRY_EUCLIDEAN_HPP
#define CLASSMAG_GEOMETRY_EUCLIDEAN_HPP

#include <array>
#include <utility>
#include <math.h> // sqrt

namespace classmag::geometry{

    template<unsigned int dimension>
    class Euclidean : public std::array<double, dimension>{
        public:
        void fill(const double d){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] = d;
        }

        double sum_(){
            auto result = 0.0;
            for (auto ii = 0u; ii < dimension; ++ii)
                result += (*this)[ii];
            return result;
        }

        Euclidean& operator=(const Euclidean &src){
            for (unsigned int ii = 0; ii < dimension; ++ii)
                (*this)[ii] = src[ii];
            return *this;
        }

        Euclidean& operator=(const std::array<double,dimension> &src){
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
        const Euclidean<dimension> &e,
        const double d){
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
    double norm(const Euclidean<dimension> &e){
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

    inline double epsilonContraction(
        unsigned int ii, 
        const Euclidean<3> &v1, 
        const Euclidean<3> &v2)
    {
        double value;
        value = v1[(ii + 1) % 3] * v2[(ii + 2) % 3] -
                v2[(ii + 1) % 3] * v1[(ii + 2) % 3];
        return value;
    };

    inline Euclidean<3> cross(
        const Euclidean<3> &e1,
        const Euclidean<3> &e2){
            Euclidean<3> w = {
                epsilonContraction(0, e1, e2),
                epsilonContraction(1, e1, e2),
                epsilonContraction(2, e1, e2)
                };
            return w;
        }

    template<unsigned int dimension>
    void circShift(Euclidean<dimension> &e){
            auto e_copy = e;
        for (unsigned int ii = 0; ii < dimension; ++ii){
            e[ii] = e_copy[(ii + 1) % dimension];
        }
    }

    template<unsigned int dimension>
    Euclidean<dimension> elementwise(
        const Euclidean<dimension> &e1, 
        const Euclidean<dimension> &e2){
        Euclidean<dimension> result;
        for (auto ii = 0u; ii < dimension; ++ii){
            result[ii] = e1[ii]*e2[ii];
        }
        return result;
    }

    template<typename T, unsigned int dimension>
    Euclidean<dimension> elementwise(
        const Euclidean<dimension> &e, 
        const std::array<T, dimension> &array){
        Euclidean<dimension> result;
        for (auto ii = 0u; ii < dimension; ++ii){
            result[ii] = e[ii]*static_cast<double>(array[ii]);
        }
        return result;
    }

    template<typename T, unsigned int dimension>
    Euclidean<dimension> elementwise(
        const std::array<T, dimension> &array,
        const Euclidean<dimension> &e){
        Euclidean<dimension> result;
        for (auto ii = 0u; ii < dimension; ++ii){
            result[ii] = e[ii]*static_cast<double>(array[ii]);
        }
        return result;
    }

    template<unsigned int dimension, unsigned int n_coeffs>
    Euclidean<dimension> linearCombination(
        const Euclidean<n_coeffs> &coefficients,
        const std::array<Euclidean<dimension>,n_coeffs> &vectors){
        Euclidean<dimension> result;
        result.fill(0.0);
        for (auto ii = 0u; ii < n_coeffs; ++ii){
            result += coefficients[ii]*vectors[ii];
        }
        return result;
    }

    template<unsigned int dim>
    bool is_shorter_than(const Euclidean<dim>& e1, const Euclidean<dim>& e2)
    {
        return (e1*e1) < (e2*e2);
    }
}



#endif