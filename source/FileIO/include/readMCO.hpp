#ifndef CLASSMAG_FILEIO_READMCO_HPP
#define CLASSMAG_FILEIO_READMCO_HPP

#include <stdexcept>
#include <map>
#include <fstream>
#include <string>
#include <iostream>

#include "Base/include/numerics.hpp"
#include "MonteCarlo/include/mcProfile.hpp"
#include "FileIO/include/stringExtractions.hpp"

namespace classmag::fileio{
    enum class MCOptions{
        THERMALIZATIONS,
        SKIPS,
        MEASUREMENTS,
        OVERRELAXATIONS,
        RESAMPLES,
        ORDERPARAMETERS
    };

    std::map<const std::string, MCOptions> mco_map_setup(){
        std::map<const std::string, MCOptions> result;
        result["Thermalizations"] = MCOptions::THERMALIZATIONS;
        result["Skips"] = MCOptions::SKIPS;
        result["Measurements"] = MCOptions::MEASUREMENTS;
        result["Overrelaxations"] = MCOptions::OVERRELAXATIONS;
        result["Resamples"] = MCOptions::RESAMPLES;
        result["Order parameters"] = MCOptions::ORDERPARAMETERS;
        return result;
    }

    int readMCO(
        montecarlo::VectorModel_Profile& target, 
        const std::string& filename
        ){
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening MCO file " + filename;
                throw std::runtime_error(errormsg); 
            }
            auto string_map = mco_map_setup();
            while (ifp.good()){
                std::string line;
                getline(ifp, line);
                auto entry = readEntryName(line);
                auto it = string_map.find(entry);
                if (it != string_map.end()){
                    switch (string_map.at(entry))
                    {
                    case MCOptions::THERMALIZATIONS:{
                        target.thermalization_ = fileio::readValue<unsigned int>(line);
                        break;
                    }
                    case MCOptions::SKIPS:{
                        target.skips_ = fileio::readValue<unsigned int>(line);
                        break;
                    }
                    case MCOptions::MEASUREMENTS:{
                        target.measurement_ = fileio::readValue<unsigned int>(line);
                        break;
                    }
                    case MCOptions::OVERRELAXATIONS:{
                        target.overrelax_ = fileio::readValue<unsigned int>(line);
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
            abort();
        }
        
        return 0;
    }

    std::vector<double> read_temperatures(const std::string& filename){
        std::vector<double> result;
        try
        {
            std::ifstream ifp(filename);
            if (!ifp.is_open() || !ifp.good()){
                std::string errormsg = "Error opening temperature file " + filename;
                throw std::runtime_error(errormsg); 
            }
            
            bool temperatures_found = false;
            while (ifp.good() && !temperatures_found){
                std::string line;
                getline(ifp, line);
                auto entry = readEntryName(line);
                if (!entry.compare("Temperatures")){
                    result = read_list<double>(line);
                    temperatures_found = true;
                }
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        return result;
    }

    #if 0
    std::vector<double> convert_temperature_vector(
        const unsigned int target_size,
        const std::vector<double>& temperatures,
        const std::vector<double>& heat_capacity
    ){

    }
    #endif
    
    void append_vector(
        std::vector<double>& target, 
        const std::vector<double>& appendix)
    {
        for (auto x : appendix)
            target.push_back(x);
    }

    std::vector<double> convert_temperature_intervals(
        const int target_size, 
        const std::vector<double>& temperature_intervals
    ){
        std::vector<double> result;
        int n_per_interval = target_size/(temperature_intervals.size() - 1);
        if (target_size > temperature_intervals.size())
        {
            int n_extra = target_size - n_per_interval * (temperature_intervals.size() - 1);
            for (auto ii = 0u; ii < temperature_intervals.size() - 1; ++ii){
                --n_extra;
                auto n_leftover = 0u;
                if (n_extra >= 0)
                    n_leftover = 1u;
                auto sub_interval = 
                    base::linspace_open(temperature_intervals[ii], temperature_intervals[ii + 1], 
                    n_per_interval + n_leftover
                    );
                append_vector(result, sub_interval);
            }
        }
        else
        {
            for (auto ii = 0u; ii < target_size; ++ii)
                result.push_back(temperature_intervals[ii]);
        }
        return result;
    }
}
#endif