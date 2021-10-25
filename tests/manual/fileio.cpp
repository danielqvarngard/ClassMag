#include <iostream>

#include "FileIO/include/optionsReadin.hpp"

using namespace classmag;

int main(int argc, char* argv[]){
    auto inputfile = fileio::get_cmd_flag_str(argc, argv, "-input");
    if (inputfile.empty())
        std::cout << "Yeehaw\n";
    return 0;
}