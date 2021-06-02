#ifndef CLASSMAG_GEOMETRY_POINTMETRIC_HPP
#define CLASSMAG_GEOMETRY_POINTMETRIC_HPP

namespace classmag::geometry{
    class PointMetric{
        public:
        virtual double distance_(const unsigned int ii, const unsigned int jj) const{
            return 0.0;
        }

        virtual double squareDistance_(
            const unsigned int site1, 
            const unsigned int site2) const {
            return 0.0;
        }

        protected:
        PointMetric(){

        }
    };
}

#endif