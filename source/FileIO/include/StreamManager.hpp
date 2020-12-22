#include <vector>
#include <string>
#include <fstream>

namespace classmag::fileio{

    class OStreamManager{
        private:
        std::ofstream os_;

        public:
        std::string newline_ = "\n";
        void newline_();
        std::string delimiter_ = ",\t";
        void delimiter_();

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
        }
    };

    class IStreamManager{
        private:
        std::ifstream is_;

        public:
        std::string delimiter_ = ",\t";
        std::string newline_ = "\n";

    };

}