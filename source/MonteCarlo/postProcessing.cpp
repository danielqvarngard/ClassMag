#include "include/postProcessing.hpp"

namespace classmag::montecarlo{
    double mean(const std::vector<double> &x){
        //Kahan summation
        auto sum = 0.0;
        auto c = 0.0;
        for (auto d : x)
        {
            auto y = d - c;
            auto t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        return sum/static_cast<double>(x.size());
    }

    double moment(const std::vector<double> &x, const double exponent){
        auto a = 0.0;
        for (auto d : x)
            a += std::pow(d, exponent);
        return a/static_cast<double>(x.size() - 1);
    }

    double variance(const std::vector<double> &x, double x_mean){
        auto sum = 0.0;
        auto c = 0.0;
        for (auto d : x)
        {
            auto x = (d - x_mean)*(d - x_mean);
            auto y = x - c;
            auto t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        sum /= static_cast<double>(x.size() - 1);
        return sum;
    }

    double variance(const std::vector<double> &x){
        auto x_mean = mean(x);
        return variance(x,x_mean);
    }

    std::pair<double, double> bootstrap(
        const std::vector<double> &data,
        std::function<double(const std::vector<double> &x)> estimator,
        unsigned int n_resamples){
            std::vector<double> estimate(n_resamples);
            
            std::mt19937 mt;
            std::uniform_int_distribution<unsigned int> distr(0, data.size() - 1);
            for (unsigned int ii = 0; ii < n_resamples; ++ii){
                std::vector<double> resample(data.size());
                for (unsigned int jj = 0; jj < data.size(); ++jj)
                    resample[jj] = data[distr(mt)];
                estimate[ii] = estimator(resample);
            }
            auto meanEstimate = mean(estimate);
            auto varianceEstimate = variance(estimate,meanEstimate);
            std::pair<double, double> result = {meanEstimate, varianceEstimate};
            return result;
    }


    DefaultBootstrapOut defaultBootstrap(
        const std::vector<double> &x, 
        const unsigned int resamples){
        auto size = x.size();

        std::uniform_int_distribution<unsigned int> distr(0,size - 1);
        std::mt19937 rng(std::random_device{}());
        std::vector<double> means(resamples);
        std::vector<double> variances(resamples);
        for (auto resample = 0u; resample < resamples; ++resample){
            std::vector<double> y(size);
            for (auto entry = 0u; entry < size; ++entry)
                y[entry] = x[distr(rng)];
            means[resample] = mean(y);
            variances[resample] = variance(y);
        }

        DefaultBootstrapOut result;
        result.meanEstimate = mean(means);
        result.meanDeviation = std::sqrt(variance(means));
        result.varianceEstimate = mean(variances);
        result.varianceDeviation = std::sqrt(variance(variances));

        return result;
    }  
}