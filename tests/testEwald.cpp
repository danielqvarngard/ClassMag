#include <iostream>

#include "Base/include/dipole.hpp"
using namespace classmag;

template<unsigned int dim>
void printEuclidean(geometry::Euclidean<dim> &e){
    for (auto d : e)
        std::cout << d << " ";
    std::cout << "\n";
}

int main(int argc, char* argv[]){
    auto swf = base::integerSweepExclude<3>(3u);
    for (auto e : swf)
        printEuclidean(e);
}