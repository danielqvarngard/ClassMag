#include <vector>
#include <string>
#include <fstream>

namespace classmag::fileio{

    class StreamManager{
        private:
        std::ofstream os_;
        std::string delimiter_ = ",\t";

        public:
        StreamManager();
        StreamManager(std::string &filename);
        StreamManager(const StreamManager &sm) = delete;
        ~StreamManager();
    };
}