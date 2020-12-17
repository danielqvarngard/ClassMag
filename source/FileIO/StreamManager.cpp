#include "include/StreamManager.hpp"

namespace classmag::fileio{
    StreamManager::StreamManager(){

    }
    StreamManager::~StreamManager(){
        os_.close();
    }
}