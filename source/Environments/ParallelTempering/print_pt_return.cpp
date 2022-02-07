#include "print_pt_return.hpp"

namespace classmag::environments
{
    void print_pt_acceptance_rate(PT_Return& source, std::ofstream& fp){
        fp << "Acceptance rates: ";
        for (auto x : source.acceptance_rates)
            fp << x << " ";
        fp << ";\n\n";
    }

    void print_pt_return_energy_only(PT_Return& source, std::ofstream& fp){
        for (auto ii = 0; ii < source.microstate_energies.size(); ++ii){
            for (auto x : source.microstate_energies[ii])
                fp << x << " ";
            fp << ";\n";
        }
    }

    void assert_variable_sizes_equal(PT_Return& source){
        if (source.is_energy_and_variables_equal_size() == false){
            std::string error_string = "PT_Return has unequal fields, ";
            std::string size_comparison = "energy size = " + source.microstate_energies.size();
            size_comparison += ", variable size = " + source.microstate_variables.size();
            error_string.append(size_comparison);
            throw std::runtime_error(error_string);
            abort();
        }
    }

    void print_pt_return(PT_Return& source, std::ofstream& fp){
        // If sizes aren't equal, we can't trust the data:
        assert_variable_sizes_equal(source);
            
        for (auto ii = 0; ii < source.microstate_energies.size(); ++ii){
            for (auto jj = 0; jj < source.microstate_energies.size(); ++jj){
                fp << source.microstate_energies[ii][jj];
                print_vector(source.microstate_variables[ii].data_[jj],fp);
                fp << "; ";
            }
            fp << "\n";
        }
    }

    void print_vector(std::vector<double>& source, std::ofstream& fp){
        for (auto x : source)
            fp << " " << x;
    }
}