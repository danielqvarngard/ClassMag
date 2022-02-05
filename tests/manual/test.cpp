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
using namespace classmag;

void print_euclidean(const geometry::Euclidean<3>& e){
    for (auto ii = 0u; ii < 3; ++ii)
        std::cout << e[ii] << " ";
    std::cout << "\n";
}

    template<unsigned int spinDim>
    void print_easyplane_norms(const montecarlo::EasyPlaneVectors<spinDim>& x){
        for (auto ii = 0; ii < x.size(); ++ii){
            std::cout << geometry::norm(x[ii][0]) << " " << geometry::norm(x[ii][1]) << "\n";
        }
    }

int main(int argc, char* argv[]){
    auto x = montecarlo::compute_tsai_easyplanes();
    print_easyplane_norms(x);
    return 0;
}