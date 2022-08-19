#ifndef CLASSMAG_PARALLELISM_ARRAYMESSAGE_HPP
#define CLASSMAG_PARALLELISM_ARRAYMESSAGE_HPP

#include <vector>

namespace classmag::parallelism{
    struct ArrayMessage{
        public:

        ArrayMessage()
        {
            
        }

        ArrayMessage(const unsigned int n):
        messageLength_(n)
        {
        }

        std::vector<std::vector<double>> data_;
        unsigned int messageLength_;
    };
}

#endif