#include "couplingLookup.hpp"
#include "Geometry/include/euclidean.hpp"

namespace classmag::base{
    template<unsigned int spinDimension> 
    class Interaction{
        public:
        virtual CouplingLookup generateTable(
            geometry::Lattice<dimension> lat){
                CouplingLookup cl();
        };
        virtual CouplingLookup generateTable(
            geometry::Lattice<dimension> lat, double cutoff);

        virtual geometry::Euclidean<spinDimension>
        protected:
        Interaction(){

        }
    };
}