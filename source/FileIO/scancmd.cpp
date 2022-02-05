#include "include/scancmd.hpp"

namespace classmag::fileio{
    int get_cmd_flag_entry(
        const int argc, 
        char* argv[], 
        const char* flag,
        const int start)
    {
        auto cont = true;
        int ii = start - 1;
        
        while (cont){
            ++ii;
            auto found = !strcmp(argv[ii],flag);
            cont = (ii < argc && !found);
        }
        return ii;
    }
        

    std::string get_cmd_flag_str(
        const int argc,
        char* argv[], 
        const char* flag)
    {
        std::string str;
        for (auto ii = 1; ii < argc; ++ii){
            if (!strcmp(argv[ii - 1], flag)){
                str = argv[ii];
            }
        }

        return str;
    }

    std::vector<double> get_cmd_vector(
        const int argc, 
        char* argv[], 
        const char* flag)
    {
        const char* delimiter_start = "[";
        const char* delimiter_end = "]";
        auto index_start = get_cmd_flag_entry(
            argc,
            argv,
            flag,
            0
        );

        if (!strcmp(argv[index_start + 1], delimiter_start)){
            // U DON GOOF
        }

        auto index_end = get_cmd_flag_entry(
            argc,
            argv,
            delimiter_end,
            index_start
        );
        std::vector<double> result;
        {
            auto n_elems = index_end - index_start - 1;
            result.resize(n_elems);
        }

        for (auto ii = index_start + 2; ii < index_end; ++ii){
            result[ii - 2] = std::stod(argv[ii]);
        }

        return result;

    }
}