#ifndef CLASSMAG_MONTECARLO_POSTPROCESSING_HPP
#define CLASSMAG_MONTECARLO_POSTPROCESSING_HPP

#include <random>
#include <vector>
#include <utility>
#include <functional>

namespace classmag::montecarlo{

    double mean(const std::vector<double> &x);
    double moment(const std::vector<double> &x, const double exponent);
    double variance(const std::vector<double> &x, double x_mean);
    double variance(const std::vector<double> &x);

    std::pair<double, double> = bootstrap(
        const std::vector<double> &data, 
        std::function<double(const std::vector<double> &)> estimator,
        const unsigned int n_resamples
        );
    
}

#endif