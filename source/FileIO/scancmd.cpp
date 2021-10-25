#include "include/scancmd.hpp"

namespace classmag::fileio{
    std::string get_cmd_flag_str(
        const int argc,
        const char* argv[], 
        const char* flag)
    {
        std::string str;
        for (auto ii = 0; ii < argc; ++ii){
            if (!strcmp(argv[ii], flag)){
                str = argv[ii + 1];
            }
        }

        return str;
    }
}