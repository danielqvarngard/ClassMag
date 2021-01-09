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

    void OStreamManager::printDelimiter_(){
        (*this) << delimiter_;
    }

    void OStreamManager::printNewline_(){
        (*this) << newline_;
    }

    std::vector<double> IStreamManager::get_CSV_row_(){
        std::vector<double> result;
        std::string line, entry;
        std::getline(is_, line);
        
        size_t pos = 0;
        while ((pos = line.find(delimiter_)) != std::string::npos){
            entry = line.substr(0,pos);
            result.push_back(std::stod(entry));
            line.erase(0, pos + delimiter_.length());
        }
        return result;
    }
}