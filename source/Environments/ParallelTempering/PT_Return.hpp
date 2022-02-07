#ifndef CLASSMAG_ENVIRONMENTS_PTRETURN_HPP
#define CLASSMAG_ENVIRONMENTS_PTRETURN_HPP

#include <vector>
#include "Parallelism/include/ArrayMessage.hpp"

namespace classmag::environments
{
    struct PT_Return
    {
        std::vector<std::vector<double>> microstate_energies;
        std::vector<parallelism::ArrayMessage> microstate_variables;
        std::vector<double> acceptance_rates;

        inline bool is_measurement_number_equal() const noexcept
        {
            return (microstate_energies.size() == microstate_variables.size());
        }

        inline bool is_process_number_equal() const noexcept
        {
            (microstate_energies[0].size() == microstate_variables[0].data_.size());
        }

        inline bool is_energy_and_variables_equal_size() const noexcept
        {
            return (
                is_measurement_number_equal() &&
                is_process_number_equal()
            );
        }
    };

}

#endif