#include "include/couplingLookup.hpp"

namespace classmag::base{
    base::CouplingLookup::CouplingLookup(
        const unsigned int n_sites, 
        const std::function<double(const int,const int)> interaction):
        n_sites_(n_sites)
        {
            couplingTable_.resize(n_sites*n_sites);
        
            for (unsigned int ii = 0; ii < n_sites*n_sites - 1; ++ii){
                int site1 = ii/n_sites;
                int site2 = ii % n_sites;
                couplingTable_[ii] = interaction(site1,site2);
            }

    }
    double CouplingLookup::coupling_(unsigned int site1, unsigned int site2) const{
        return couplingTable_[n_sites_*site1 + site2];
    }
}