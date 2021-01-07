#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cerrno>

namespace classmag::fileio{

    class OStreamManager{
        private:
        std::ofstream os_;

        public:
        std::string delimiter_ = ",\t";
        void delimiter_();
        std::string newline_ = "\n";
        void newline_();

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
    };

    class IStreamManager{
        private:
        std::ifstream is_;

        public:
        std::string delimiter_ = ",\t";
        std::string newline_ = "\n";

        IStreamManager &operator>>(double &x);

        std::vector<double> get_CSV_row_();

        template<unsigned int size>
        IStreamManager &operator>>(std::array<double,size> &array){
            std::string headers;
            std::getline(is_,headers);
            std::stringstream entryBuffer(headers);
            std::string entry;
            try{
                for (auto ii = 0; ii < array.size(); ++ii){
                    (*this) >> array[ii];
                }
            }
        }
    };

}