#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <mpi.h>

int main(int argc, char* argv[]){
    std::ifstream ifp;
    for (auto ii = 0; ii < argc; ++ii){
        if (!std::strcmp(argv[ii], "-input")){
            std::cout << "Input flag detected\n";
            ifp.open(argv[ii + 1], std::ifstream::in);
        }
    }

    int s = 0;
    int r = 0;

    if (ifp.is_open()){
        std::cout << "ifp.is_open()\n";
        
        while (ifp.good()){
            std::string str;
            std::stringstream str_stream;
            getline(ifp, str);
            std::cout << str << "\n";
            str_stream << str;
            while (!str_stream.eof()){
                std::string temp_str;
                str_stream >> temp_str;
                int temp_int;
                if (std::stringstream(temp_str) >> r)
                    std::cout << "r = " << r << std::endl;
            }

        }
    }
    return 0;
}