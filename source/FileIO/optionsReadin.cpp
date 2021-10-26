#include "include/optionsReadin.hpp"

namespace classmag::fileio{

    std::string get_input_file_name(int argc, char* argv[]){
        auto filename = get_cmd_flag_str(argc, argv, "-input");
        return filename;
    }

    std::ifstream get_input_stream(int argc, char* argv[]){
        auto filename = get_input_file_name(argc, argv);
        std::ifstream ifp;
        try
            {
                ifp.open(filename, std::ifstream::in);
                if(!ifp.is_open() || !ifp.good())
                    throw std::runtime_error("Input file not found.");
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << std::endl;
                abort();
            }
        return ifp;
    }
    
    
}
