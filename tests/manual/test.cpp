#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <array>
#include <vector>
#include <map>
#include <stdexcept>

#include "Geometry/include/predefLattices.hpp"
#include "MonteCarlo/include/PredefAnisotropyVectors.hpp"
#include "Base/include/numerics.hpp"
#include "FileIO/include/readMCO.hpp"
using namespace classmag;

template<typename T>
geometry::Euclidean<3> convert_array(const std::array<T, 3>& a)
{
    geometry::Euclidean<3> e;
    for (auto ii = 0u; ii < 3; ++ii)
    {
        e[ii] = static_cast<double>(a[ii]);
    }
    return e;
}

void print_euclidean(const geometry::Euclidean<3>& e){
    for (auto ii = 0u; ii < 3; ++ii)
        std::cout << e[ii] << " ";
    std::cout << "\n";
}

void print_euclidean(const std::vector<geometry::Euclidean<3>>& v){
    for (auto s : v)
        print_euclidean(s);
}

template<unsigned int spinDim>
void print_easyplane_norms(const montecarlo::EasyPlaneVectors<spinDim>& x){
    for (auto ii = 0; ii < x.size(); ++ii){
        std::cout << geometry::norm(x[ii][0]) << " " << geometry::norm(x[ii][1]) << "\n";
    }
}



int main(int argc, char* argv[]){
    for (auto ii = 1u; ii < base::maximum_z3plus_index(3); ++ii)
    {
        auto a = base::z3plus_entry(ii,3);
        auto e = convert_array(a);
        print_euclidean(e);
    }

    return 0;
}