#ifndef CLASSMAG_FILEIO_FILESYSTEM_HPP
#define CLASSMAG_FILEIO_FILESYSTEM_HPP

#include <ctime> 
#include <string>

namespace classmag::fileio{
    std::string datestamp();

    struct VectorTarget{
        public:
        VectorTarget(const unsigned int n):
        messageLength_(n)
        {

        }

        std::vector<std::vector<double>> data_;
        const unsigned int messageLength_;
    };
}

#endif