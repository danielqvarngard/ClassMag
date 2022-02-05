#ifndef CLASSMAG_FILEIO_FILESYSTEM_HPP
#define CLASSMAG_FILEIO_FILESYSTEM_HPP

#include <ctime> 
#include <string>
#include <vector>
#include <fstream>

namespace classmag::fileio{
    std::string datestamp();


    struct DefaultMonteCarloParameters
    {
        int thermalizations;
        int skipped_sweeps;
        int measurements;
        int overrelaxations;
        int order_parameters;

        double nearest_neighbor_strength;
        double nearest_neighbor_cutoff;
        
        double rkky_strength;
        double rkky_wavevector;
        double rkky_cutoff;

        double dipole_strength;
        double dipole_k_real;
        double dipole_k_reciprocal;
        double dipole_alpha = -1.0;

        std::vector<int> system_size;
    };

    DefaultMonteCarloParameters ScanDefaultMontecarloParameters(
        const std::string infile
    );

    void WriteDefaultFileHeader(
        const std::string& outfile, 
        const DefaultMonteCarloParameters& parameters
    );
}

#endif