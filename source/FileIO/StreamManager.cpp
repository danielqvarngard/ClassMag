#include "include/StreamManager.hpp"

namespace classmag::fileio{
    OStreamManager::OStreamManager():
    os_(std::ofstream("filestream.out"))
    {

    }

    OStreamManager::OStreamManager(const std::string &filename):
    os_(std::ofstream(filename))
    {

    }

    OStreamManager::~OStreamManager(){
        os_.close();
    }

    void OStreamManager::delimiter_(){
        (*this) << delimiter_;
    }

    void OStreamManager::newline_(){
        (*this) << newline_;
    }
}