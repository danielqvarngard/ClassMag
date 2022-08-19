#include "include/ParseOutputValues.hpp"

namespace classmag::fileio{
    void push_line(NaiveSimulationDataContainer& collection, const std::string& line)
    {
        
        auto streamer = std::stringstream(line);
        while (!streamer.eof() && streamer.good()){
            std::vector<double> numbers;
            double holder;
            std::string variable_set;
            std::getline(streamer, variable_set, ';'); 
            std::stringstream variable_streamer(variable_set);
            while (variable_streamer >> holder)
                numbers.push_back(holder);
        
            for (auto ii = 0; ii < numbers.size(); ++ii){
                collection.variables[ii].push_back(numbers[ii]);
            }
        }
    }

    void read_results(NaiveSimulationDataContainer& collection, std::ifstream& file_pointer)
    {
        bool there_is_more_data = true;
        while (there_is_more_data)
        {
            std::string line;
            std::getline(file_pointer, line);
            bool is_not_at_matrix_end = (line.find('}') == std::string::npos);

            if (is_not_at_matrix_end)
            {
                push_line(collection,line);
            }

            there_is_more_data = 
                is_not_at_matrix_end &&
                file_pointer.good() &&
                !file_pointer.eof();
        }
    }

    std::ifstream& scan_for_start_of_results(std::ifstream& ifp)
    {
        bool results_found = false;
        while (ifp.good() && !results_found){
            std::string line;
            getline(ifp, line);
            auto entry = readEntryName(line);
            if (!entry.compare("Results")){
                results_found = true;
            }
        }

        return ifp;
    }
}