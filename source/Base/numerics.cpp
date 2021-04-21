#include "include/numerics.hpp"

namespace classmag::base{
    double pi(){
        return M_PI;
    }

    std::vector<double> shanks(std::vector<double> &x){
        std::vector<double> s(x.size() - 2);
        for (auto ii = 0u; ii < s.size(); ++ii){
            auto num = (x[ii + 2] - x[ii + 1]) * (x[ii + 2] - x[ii + 1]);
            auto denom = x[ii + 2] - 2.0 * x[ii + 1] + x[ii];
            s[ii] = x[ii + 2] - num/denom;
        }
        return s;
    }
}