#include "include/filesystem.hpp"

namespace classmag::fileio{
    std::string datestamp(){
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,sizeof(buffer),"%Y%m%d_%H-%M",timeinfo);
        std::string stamp(buffer);
        return stamp;
    }
}