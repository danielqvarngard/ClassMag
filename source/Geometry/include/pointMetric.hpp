#ifndef CLASSMAG_GEOMETRY_POINTMETRIC_HPP
#define CLASSMAG_GEOMETRY_POINTMETRIC_HPP

namespace classmag::geometry{
    class PointMetric{
        public:
        virtual double distance(const unsigned int ii, const unsigned int jj) const{
            return 0.0;
        }

        protected:
        PointMetric(){

        }
    };
}

#endif