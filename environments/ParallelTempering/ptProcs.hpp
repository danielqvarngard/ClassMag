#ifndef CLASSMAG_ENVIRONMENTS_PTPROCS_HPP
#define CLASSMAG_ENVIRONMENTS_PTPROCS_HPP

#include "MonteCarlo/include/VectorModelManager.hpp"
namespace classmag::environments{
    int mpiPT_hub(unsigned int seed);
    int mpiPT_mc();
}

#endif