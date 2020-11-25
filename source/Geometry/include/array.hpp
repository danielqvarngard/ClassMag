#ifndef CLASSMAG_GEOMETRYLIB_ARRAY_HPP
#define CLASSMAG_GEOMETRYLIB_ARRAY_HPP

#include <vector>
#include <math.h>

namespace classmag::geometry{

    class Array{
        private:
            std::vector<double> elements_;
        public:
            const unsigned int rows_;
            const unsigned int columns_;
            Array(
                const unsigned int rows, 
                const unsigned int columns);
            Array(const Array &array);
            ~Array();
            Array(const Array &&array);
            void set(const std::vector<double> target);

            
    };
}

#endif