#ifndef CLASSMAG_PARALLELISM_ARRAYMESSAGE_HPP
#define CLASSMAG_PARALLELISM_ARRAYMESSAGE_HPP

#include <vector>

namespace classmag::parallelism{
    struct VectorTarget{
        public:
        VectorTarget(const unsigned int n):
        messageLength_(n)
        {
            data_.resize(n);
        }

        std::vector<std::vector<double>> data_;
        const unsigned int messageLength_;
    };
}

#endif