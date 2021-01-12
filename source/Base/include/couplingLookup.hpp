#ifndef CLASSMAG_MODELS_COUPLINGLOOKUP_HPP
#define CLASSMAG_MODELS_COUPLINGLOOKUP_HPP
#include <vector>
#include <functional>
#include <iostream>

namespace classmag::models{
    class CouplingLookup{
        private:
            std::vector<double> couplingTable_;
            const unsigned int n_sites_;
        public:
            CouplingLookup(
                const unsigned int n_sites, 
                const std::function<double(int, int)> interaction);
            double coupling_(unsigned int site1, unsigned int site2) const;
    };
}

#endif