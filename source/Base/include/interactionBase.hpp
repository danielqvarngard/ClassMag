#include "couplingLookup.hpp"
#include "Geometry/include/lattice.hpp"

namespace classmag::base{
    template<unsigned int dimension>
    class InteractionBase{
        public:
        virtual CouplingLookup generateTable(
            geometry::Lattice<dimension> lat){
                CouplingLookup cl();
        };
        virtual CouplingLookup generateTable(
            geometry::Lattice<dimension> lat, double cutoff);
        protected:
        InteractionBase(){

        }
    };
}