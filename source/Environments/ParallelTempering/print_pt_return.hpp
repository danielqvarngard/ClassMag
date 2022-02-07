#ifndef CLASSMAG_ENVIRONMENTS_PRINT_PT_RETURN_HPP
#define CLASSMAG_ENVIRONMENTS_PRINT_PT_RETURN_HPP

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>

#include "PT_Return.hpp"

namespace classmag::environments
{
    void print_pt_acceptance_rate(PT_Return& source, std::ofstream& fp);
    void print_pt_return_energy_only(PT_Return& source, std::ofstream& fp);
    void assert_variable_sizes_equal(PT_Return& source);
    void print_vector(std::vector<double>& source, std::ofstream& fp);
    void print_pt_return(PT_Return& source, std::ofstream& fp);
}

#endif