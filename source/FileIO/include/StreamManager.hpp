#ifndef CLASSMAG_FILEIO_STREAMMANAGER_HPP
#define CLASSMAG_FILEIO_STREAMMANAGER_HPP

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cerrno>

#include "filesystem.hpp"
#include "Parallelism/include/ArrayMessage.hpp"

namespace classmag::fileio{
    class OStreamManager{
        private:
        std::ofstream os_;

        public:
        std::string delimiter_ = ",\t";
        void printDelimiter_();
        std::string newline_ = "\n";
        void printNewline_();

        OStreamManager();
        OStreamManager(const std::string &filename);
        OStreamManager(const OStreamManager &sm) = delete;
        ~OStreamManager();

        template <typename T>
        OStreamManager &operator<<(const T x){
            os_ << x;
            return (*this);
        }

        template <typename T>
        OStreamManager &operator<<(const std::vector<T> &v){
            for (auto x : v){
                (*this) << x;
                (*this) << delimiter_;
            }
            (*this) << newline_;        
            return (*this);
        }

        template <typename T, unsigned int size>
        OStreamManager &operator<<(const std::array<T,size> &v){
            for (auto x : v){
                (*this) << x;
                (*this) << delimiter_;
            }
            (*this) << newline_;
            return (*this);
        }

        OStreamManager &operator<<(const parallelism::ArrayMessage &v){
            for (auto ii = 0u; ii < v.messageLength_; ++ii){
                for (auto d : v.data_[ii]){
                    (*this) << d;
                    (*this) << delimiter_;
                }
            }
            (*this) << newline_;
            return (*this);
        }
    };
}

#endif