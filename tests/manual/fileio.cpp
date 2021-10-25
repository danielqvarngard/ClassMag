#include <iostream>

#include "FileIO/include/optionsReadin.hpp"

using namespace classmag;

int main(int argc, char* argv[]){
    auto file_name = fileio::get_cmd_flag_str(argc, argv, "-input");
    if (file_name.empty())
        std::cout << "Yeehaw\n";
    return 0;
}