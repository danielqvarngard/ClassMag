#ifndef CLASSMAG_FILEIO_STRINGEXTRACTIONS
#define CLASSMAG_FILEIO_STRINGEXTRACTIONS

#include <string>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <stdexcept>
#include <map>

namespace classmag::fileio{

    inline std::string breakstring(){
        return "};";
    }

    std::string readEntryName(const std::string &row){
        std::stringstream stream;
        stream << row;
        bool finito = false;
        std::string entry;
        while (!stream.eof() && !finito){
            
            std::string temp;
            stream >> temp;
            if (temp.compare("=") != 0){
                if (!entry.empty())
                    entry += " ";
                entry += temp;
            }
            else
                finito = true;
        }
        return entry;
    }

    int readDelimiter(const std::string &line){
        if (line.find('{') < line.npos)
            return -1;
        if (line.find('}') < line.npos)
            return +1;
        return 0;
    }

    bool matrixEnd(const std::string &line){
        if (line.find('}') < line.npos)
            return true;
        else
            return false;
    }

    template<typename T, unsigned int dim>
    std::array<T, dim> readRow(const std::string line){
        std::stringstream stream;
        std::array<T, dim> result;
        stream << line;
        unsigned int ii = 0;
        while (!stream.eof()){
            std::string temp_str;
            T value;
            stream >> temp_str;
            if (std::stringstream(temp_str) >> value){
                result[ii] = value;
                ++ii;
            }
        }
        if (ii < dim){
            std::cerr << "Dimension mismatch\n";
            abort(); 
        }
        return result;
        
    }

    template<typename T, unsigned int columns>
    std::vector<std::array<T, columns>> readMatrix(std::ifstream &ifp){
        std::vector<std::array<T,columns>> result;
        bool endOfMatrix = false;
        while (!endOfMatrix){
            std::string line;
            getline(ifp, line);
            if (matrixEnd(line) || ifp.eof())
                endOfMatrix = true;
            else{
                auto row = readRow<T,columns>(line);
                result.push_back(row);
            }
        }
        return result;
    }
}

#endif