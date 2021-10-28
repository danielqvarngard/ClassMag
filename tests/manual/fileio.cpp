#include <iostream>

#include "FileIO/include/optionsReadin.hpp"
#include "FileIO/include/readCouplings.hpp"

using namespace classmag;

int main(int argc, char* argv[]){
    auto file_name = fileio::get_input_file_name(argc, argv);
    base::CouplingsMatrixDense<3u> s(0);
    auto resultflag = fileio::readLinearInteractions<geometry::Matrix<3, 3>,3,3>(s, file_name);
    return 0;
}