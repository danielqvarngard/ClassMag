#include "include/postProcessing.hpp"

namespace classmag::montecarlo{
    double mean(const std::vector<double> &x){
        auto a = 0.0;
        for (auto d : x)
            a += d;
        return a/static_cast<double>(x.size());
    }

    double moment(const std::vector<double> &x, const double exponent){
        auto a = 0.0;
        for (auto d : x)
            a += std::pow(d, exponent);
        return a/static_cast<double>(x.size());
    }

    double variance(const std::vector<double> &x, double x_mean){
        auto a = moment(x,2.0);
        return a - x_mean*x_mean;
    }

    double variance(const std::vector<double> &x){
        auto x_mean = mean(x);
        return variance(x,x_mean);
    }

    std::pair<double, double> bootstrap(
        const std::vector<double> &data,
        std::function<double(std::vector<double> &)> estimator,
        unsigned int n_resamples){
            std::vector<double> estimate(n_resamples);
            
            std::mt19937 mt;
            std::uniform_int_distribution<unsigned int> distr(0, data.size() - 1);
            for (unsigned int ii = 0; ii < n_resamples; ++ii){
                std::vector<double> resample(data.size());
                for (unsigned int jj = 0; jj < data.size(); ++ii)
                    resample[jj] = data[distr(mt)];
                estimate[ii] = estimator(resample);
            }
            auto meanEstimate = mean(estimate);
            auto varianceEstimate = variance(estimate,meanEstimate);
            std::pair<double, double> result = {meanEstimate, varianceEstimate};
            return result;
        }
}