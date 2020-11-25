#ifndef CLASSMAG_GEOMETRY_QARTESIAN_HPP
#define CLASSMAG_GEOMETRY_QARTESIAN_HPP

#include <memory>
#include <cmath>
#include <vector>
#include <initializer_list>
#include <iostream>

namespace classmag::geometry
{
    template<unsigned int size>
    class Qartesian
    {
        private:
        std::unique_ptr<double[]> entry;

        public:
        Qartesian<size>():
        entry(std::make_unique<double[]>(size))
        {
            
        }

        Qartesian<size>(const double initValue):
        entry(std::make_unique<double[]>(size))
        {
            for (unsigned int index = 0; index < size; ++index)
                entry[index] = initValue;
        }

        Qartesian<size>(std::initializer_list<double> l):
        entry(std::make_unique<double[]>(size))
        {
            if (l.size() != size)
            {   
                throw "Qartesian initialization list error!";
            }
            unsigned int ii = 0;
            for (auto element : l)
            {
                entry[ii] = element;
                ++ii;
            }
        }

        Qartesian<size>(const Qartesian &cartesian):
        entry(std::make_unique<double[]>(size))
        {
            set(cartesian);
        }

        double get(unsigned int index) const;
        void set(unsigned int index, const double value);
        void set(const Qartesian<size> &cartesian);
        void print();
        void operator=(const Qartesian<size> &target);
    };

    template <unsigned int n>
    double Qartesian<n>::get(unsigned int index) const
    {
        double result;
        result = entry[index];
        return result;
    }

    template <unsigned int n>
    void Qartesian<n>::set(unsigned int index, const double value)
    {
        entry[index] = value;
    }

    template <unsigned int n>
    void Qartesian<n>::set(const Qartesian<n> &cartesian)
    {
        for (unsigned int ii = 0; ii < n; ++ii)
            set(ii, cartesian.get(ii));
    }
    
    template <unsigned int n>
    void Qartesian<n>::print()
    {
        for (unsigned int ii = 0; ii < n; ++ii)
            std::cout << get(ii) << "\n";
    }

    template <unsigned int n>
    void Qartesian<n>::operator=(const Qartesian<n> &target)
    {
        for (unsigned int ii = 0; ii < n; ++ii)
            entry[ii] = target.get(ii);
    }

    template <unsigned int n>
    Qartesian<n> operator*(const double magnitude, const Qartesian<n> &cartesian)
    {
        Qartesian<n> result;
        for (unsigned int ii = 0; ii < n; ++ii)
            result.set(ii, magnitude*cartesian.get(ii));
        return result;
    }

    template <unsigned int n>
    Qartesian<n> operator*(const int magnitude, const Qartesian<n> &cartesian)
    {
        Qartesian<n> result;
        for (unsigned int ii = 0; ii < n; ++ii)
            result.set(ii, static_cast<double>(magnitude)*cartesian.get(ii));
        return result;
    }

    template <unsigned int n>
    Qartesian<n> operator*(const Qartesian<n> &cartesian, const double magnitude)
    {
        return magnitude * cartesian;
    }

    template <unsigned int n>
    Qartesian<n> operator/(const Qartesian<n> &cartesian, const double magnitude)
    {
        return 1.0/magnitude * cartesian;
    }

    template <unsigned int n>
    Qartesian<n>operator+(const Qartesian<n>&cartesian1, const Qartesian<n>&cartesian2)
    {
        Qartesian<n>result;
        for (unsigned int ii = 0; ii < n; ++ii)
            result.set(ii, cartesian1.get(ii) + cartesian2.get(ii));
        return result;
    }

    template <unsigned int n>
    Qartesian<n>operator-(const Qartesian<n>&cartesian1, const Qartesian<n>&cartesian2)
    {
        Qartesian<n> result = cartesian1 + (-1.0) * cartesian2;
        return result;
    }

    template <unsigned int n>
    void operator+=(Qartesian<n> &cartesian1, const Qartesian<n> &cartesian2)
    {
        cartesian1 = cartesian1 + cartesian2;
    }

    template <unsigned int n>
    double operator*(const Qartesian<n> &cartesian1, const Qartesian<n> &cartesian2)
    {
        double result = 0.0;
        for (unsigned int ii = 0; ii < n; ++ii)
            result += cartesian1.get(ii)*cartesian2.get(ii);
        return result;
    }


    template <unsigned int n>
    double norm(const Qartesian<n> &cartesian)
    {
        double norm = std::sqrt((cartesian) * (cartesian));
        return norm;
    }

    template <unsigned int n>
    void normalize(Qartesian<n> &cartesian)
    {
        cartesian = cartesian/norm(cartesian);
    }

    template <unsigned int n>
    Qartesian<n> proj(const Qartesian<n> &vector, const Qartesian<n> &direction)
    {
        Qartesian<n> result;
        try
        {
            result = (vector * direction)/(direction * direction) * direction;
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        return result;
    }

    template <unsigned int n>
    Qartesian<n> projNormalized(const Qartesian<n> &vector, const Qartesian<n> &unitDirection)
    {
        Qartesian<n> result;
        result = (vector * unitDirection) * unitDirection;
        return result;
    }

    inline double epsilonContraction(unsigned int ii, const Qartesian<3> &v1, const Qartesian<3> &v2)
    {
        double value;
        value = v1.get((ii + 1) % 3) * v2.get((ii + 2) % 3) -
                v2.get((ii + 1) % 3) * v1.get((ii + 2) % 3);
        return value;
    };

    inline Qartesian<3> cross(const Qartesian<3> &v1, const Qartesian<3> &v2)
    {
        Qartesian<3> w;
        for (unsigned int ii = 0; ii < 3; ++ii)
        {   
            double value = epsilonContraction(ii,v1,v2);
            w.set(ii,value);
        }
        return w;
    }

    inline void mapToGlobal(Qartesian<3> &cartesian, const Qartesian<3> &old_z_axis)
    {
        Qartesian<3> v1 = {old_z_axis.get(2), 0.0, -old_z_axis.get(0)};
        normalize(v1);

        Qartesian<3> v2 = cross(old_z_axis, v1);

        cartesian = (cartesian.get(0) * v1) + (cartesian.get(1) * v2) + (cartesian.get(2) * old_z_axis);
    }
} // namespace classmag::geometry


#endif
