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

    std::vector<double> linspace(
        const double min, 
        const double max, 
        const unsigned int stepCount){
        
        std::vector<double> result(stepCount);
        if (stepCount < 2){
            result.resize(1);
            result[0] = min;
            return result;
        }

        auto n = static_cast<double>(stepCount);
        auto delta = (max - min)/(n - 1.0);
        for (auto ii = 0u; ii < stepCount; ++ii){
            auto x = static_cast<double>(ii);
            result[ii] = min + delta * x;
        }
        return result;
    }
}