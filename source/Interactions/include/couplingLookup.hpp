#ifndef CLASSMAG_INTERACTIONS_COUPLINGLOOKUP_HPP
#define CLASSMAG_INTERACTIONS_COUPLINGLOOKUP_HPP
#include <vector>
#include <functional>

namespace classmag::interactions{
    class CouplingLookup{
        private:
            std::vector<double> couplingTable_;
            std::vector<unsigned int> couplingColumn_;
        public:
            CouplingLookup(
                const unsigned int n_sites, 
                const std::function<double(int, int)> interaction);
            double coupling_(unsigned int site1, unsigned int site2);
    };
}

#endif