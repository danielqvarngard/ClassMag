#ifndef CLASSMAG_MONTECARLO_POSTPROCESSING_HPP
#define CLASSMAG_MONTECARLO_POSTPROCESSING_HPP

#include <random>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include "Base/include/magnetizationData.hpp"

namespace classmag::montecarlo{

    double mean(const std::vector<double> &x);
    double moment(const std::vector<double> &x, const double exponent);
    double variance(const std::vector<double> &x, double x_mean);
    double variance(const std::vector<double> &x);

    std::pair<double, double> bootstrap(
        const std::vector<double> &data,
        std::function<double(const std::vector<double> &x)> estimator,
        const unsigned int n_resamples);

    struct DefaultBoostrapOut
    {
    public:
        double meanEstimate;
        double meanDeviation;
        double varianceEstimate;
        double varianceDeviation;
    };

    template <unsigned int spinDimension, unsigned int subLattices>
    struct VectorStatEstimate{
        base::MagnetizationData<spinDimension, subLattices> mean_;
        base::MagnetizationData<spinDimension, subLattices> variance_;
    };

    template <unsigned int spinDimension, unsigned int subLattices>
    struct VectorBootstrapOut
    {
        public:
        base::MagnetizationData<spinDimension, subLattices> meanEstimate_;
        base::MagnetizationData<spinDimension, subLattices> meanVariance_;
        base::MagnetizationData<spinDimension, subLattices> varianceEstimate_;
        base::MagnetizationData<spinDimension, subLattices> varianceVariance_;
    };

    DefaultBoostrapOut defaultBootstrap(const std::vector<double> &x);    

    template <unsigned int spinDimension, unsigned int subLattices>
    base::MagnetizationData<spinDimension, subLattices>
    vectorMean(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &x){
        
        base::MagnetizationData<spinDimension, subLattices> result;
        result.fill_(0.0);

        for (auto v : x)
            result += v;

        result *= 1.0/static_cast<double>(x.size());
        return result;
    }

    template <unsigned int spinDimension, unsigned int subLattices>
    base::MagnetizationData<spinDimension, subLattices>
    vectorVariance(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &x,
        base::MagnetizationData<spinDimension, subLattices> &mean){
        
        base::MagnetizationData<spinDimension, subLattices> result;
        result.fill_(0.0);

        for (auto v : x)
            result += (v - mean) * (v - mean);

        result *= 1.0/static_cast<double>(x.size() - 1);
        return result;
    }

    template <unsigned int spinDimension, unsigned int subLattices>
    base::MagnetizationData<spinDimension, subLattices>
    vectorVariance(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &x){
            
        auto mean = vectorMean(x);
        return vectorVariance(x,mean);
    }

    template <unsigned int spinDimension, unsigned int subLattices>
    VectorBootstrapOut<spinDimension, subLattices> VectorBootstrap(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &data,
        const unsigned int resamples)
    {
        std::mt19937 mt;
        std::uniform_int_distribution<unsigned int> distr(0, data.size() - 1);

        std::vector<base::MagnetizationData<spinDimension, subLattices>> 
            meanEstimates(resamples), varianceEstimates(resamples);
        for (auto ii = 0u; ii < resamples; ++ii){
            const std::vector<base::MagnetizationData<spinDimension, subLattices>&> 
                resampleAdresses(data.size());
            for (auto jj = 0u; jj < data.size(); ++jj){
                auto randEntry = distr(mt);
                resampleAdresses[jj] = &(data[randEntry]);
            }

            meanEstimates[ii] = vectorMean(resampleAdresses);
            varianceEstimates[ii] = vectorVariance(resampleAdresses);
        }
        VectorBootstrapOut<spinDimension, subLattices> result;
        result.meanEstimate_ = vectorMean(meanEstimates);
        result.meanVariance_ = vectorVariance(meanEstimates);
        result.varianceEstimate_ = vectorMean(varianceEstimates);
        result.varianceVariance_ = vectorVariance(varianceEstimates);
        return result;
    }

    template <unsigned int spinDimension, unsigned int subLattices>
    geometry::Euclidean<subLattices> trChi(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &x, 
        const double beta
    ){
        auto variances = vectorVariance(x);
        geometry::Euclidean<subLattices> result;
        for (auto ii = 0u; ii < subLattices; ++ii){
            result[ii] = beta*variances[ii].sum_();
        }
        return result;
    }

    template <unsigned int spinDimension, unsigned int subLattices>
    std::array<std::pair<double, double>, subLattices> trChiBootstrap(
        const std::vector<base::MagnetizationData<spinDimension, subLattices>> &data,
        const double beta,
        const unsigned int resamples)
    {
        std::mt19937 mt;
        std::uniform_int_distribution<unsigned int> distr(0, data.size() - 1);

        std::array<std::vector<double>, subLattices> estimates;
        for (auto ii = 0u; ii < subLattices; ++ii)
            estimates[ii].resize(resamples);
        
        for (auto ii = 0u; ii < resamples; ++ii){
            std::vector<base::MagnetizationData<spinDimension, subLattices>> 
                resampleAdresses(data.size());
            for (auto jj = 0u; jj < data.size(); ++jj){
                auto randEntry = distr(mt);
                resampleAdresses[jj] = data[randEntry];
            }
            auto estimate = trChi(resampleAdresses, beta);
            for (auto jj = 0u; jj < subLattices; ++jj)
                estimates[jj][ii] = estimate[jj];
        }
        std::array<std::pair<double, double>, subLattices> result;
        for (auto ii = 0u; ii < subLattices; ++ii){
            result[ii].first = mean(estimates[ii]);
            result[ii].second = sqrt(variance(estimates[ii]));
        }

        
        return result;

    }
}

#endif