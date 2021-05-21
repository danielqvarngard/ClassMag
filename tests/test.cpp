#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <array>
#include <vector>

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
    if (ii < dim)
        std::cout << "error reading entries\n";
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

template<typename T, unsigned int dim>
void printrow(std::array<T, dim> &row){
    for(auto ii = 0u; ii < dim; ++ii)
        std::cout << row[ii] << " ";
}

template<typename T, unsigned int dim>
void printmatrix(std::vector<std::array<T, dim>> &matrix){
    for(auto ii = 0u; ii < matrix.size(); ++ii){
        printrow<T, dim>(matrix[ii]);
        std::cout << "\n";
    }
}

int main(int argc, char* argv[]){
    std::ifstream ifp;
    for (auto ii = 0; ii < argc; ++ii){
        if (!std::strcmp(argv[ii], "-input")){
            std::cout << "Input flag detected\n";
            ifp.open(argv[ii + 1], std::ifstream::in);
        }
    }

    int s = 0;
    int r = 0;

    if (ifp.is_open() && ifp.good()){
        std::cout << "ifp.is_open()\n";
        
        while (ifp.good()){
            std::string str;
            std::stringstream str_stream;
            getline(ifp, str);
            str_stream << str;
            std::string entry;
            if(readEntryName(str).compare("Array") == 0){
                auto arr = readMatrix<double,3>(ifp);
                printmatrix<double,3>(arr);
            }

            while (!str_stream.eof()){
                std::string temp_str;
                str_stream >> temp_str;
                double temp_value;
            }

        }
    }
    return 0;
}