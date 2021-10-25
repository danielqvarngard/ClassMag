#ifndef CLASSMAG_FILEIO_SCANCMD_HPP
#define CLASSMAG_FILEIO_SCANCMD_HPP
#include <string>
#include <sstream>

namespace classmag::fileio{
    std::string get_cmd_flag_str(
        const int argc, 
        char* argv[], 
        const char* flag);
}

#endif