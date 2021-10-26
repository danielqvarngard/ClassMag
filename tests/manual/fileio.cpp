#include <iostream>

#include "FileIO/include/optionsReadin.hpp"
#include "FileIO/include/readCouplings.hpp"

using namespace classmag;

int main(int argc, char* argv[]){
    auto file_name = fileio::get_input_file_name(argc, argv);
    base::CouplingsMatrixDense<3> s(0);
    //auto couplings = fileio::readLinearInteractions(s, file_name);
    return 0;
}