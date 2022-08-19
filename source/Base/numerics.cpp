#include "include/numerics.hpp"

namespace classmag::base{

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

    std::vector<double> linspace_open(
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
        auto delta = (max - min)/n;
        for (auto ii = 0u; ii < stepCount; ++ii){
            auto x = static_cast<double>(ii);
            result[ii] = min + delta * x;
        }
        return result;
    }

    std::vector<double> logspace(
        const double min,
        const double max,
        const unsigned int step_count
    ){
        std::vector<double> result(step_count);
        const auto c = std::pow(max/min,1/(step_count - 1));
        result[0] = min;

        for (auto ii = 1; ii < step_count; ++ii){
            result[ii] = c*result[ii - 1];
        }
        return result;
    }
}