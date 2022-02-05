#ifndef CLASSMAG_FILEIO_SCANCMD_HPP
#define CLASSMAG_FILEIO_SCANCMD_HPP
#include <string.h>
#include <sstream>
#include <vector>

namespace classmag::fileio{
    int get_cmd_flag_entry(
        const int argc, 
        char* argv[], 
        const char* flag,
        const unsigned int start);

    std::string get_cmd_flag_str(
        const int argc, 
        char* argv[], 
        const char* flag);

    std::vector<double> get_cmd_vector(
        const int argc, 
        char* argv[], 
        const char* flag);
}

#endif