#ifndef CLASSMAG_MONTECARLO_MCPROFILE_HPP
#define CLASSMAG_MONTECARLO_MCPROFILE_HPP

namespace classmag::montecarlo{
    struct MC_Profile{
        public:
        unsigned int thermalization_;
        unsigned int measurement_;
        unsigned int skips_;
        unsigned int overrelax_;
        unsigned int n_sites_;
        unsigned int seed_;
    };
}

#endif 