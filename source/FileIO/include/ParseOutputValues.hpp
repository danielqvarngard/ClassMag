#ifndef CLASSMAG_FILEIO_PARSE_OUTPUT_VALUES_HPP
#define CLASSMAG_FILEIO_PARSE_OUTPUT_VALUES_HPP

#include "stringExtractions.hpp"

#include <memory>
#include <vector>
#include <fstream>

namespace classmag::fileio{
    struct NaiveSimulationDataContainer{
        std::vector<double> temperatures;
        std::vector<std::vector<double>> variables;
    };

    void push_line(NaiveSimulationDataContainer& collection, const std::string& line);
    void read_results(NaiveSimulationDataContainer& collection, std::ifstream& ifp);
    std::ifstream& scan_for_start_of_results(std::ifstream& ifp);
    
}

#endif