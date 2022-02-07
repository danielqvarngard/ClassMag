#include "include/filesystem.hpp"

namespace classmag::fileio{
    std::string datestamp(){
        time_t rawtime;
        struct tm * timeinfo;
        char buffer[80];
        time (&rawtime);
        timeinfo = localtime(&rawtime);
        strftime(buffer,sizeof(buffer),"%Y%m%d_%H-%M",timeinfo);
        std::string stamp(buffer);
        return stamp;
    }

    void WriteDefaultFileHeader(
        std::ofstream& of,
        const DefaultMonteCarloParameters& parameters
    ){

        of << "--- Update parameters: ---\n\n";
        of << "Thermalizations = " << parameters.thermalizations << ";\n";
        of << "Skips = " << parameters.skipped_sweeps << ";\n";
        of << "Measurements = " << parameters.measurements << ";\n";
        of << "Overrelaxations = " << parameters.overrelaxations << ";\n";
        
        of << "\n--- Coupling descriptions: ---\n\n";
        of << "Nearest neighbor = {\n";
        of << "\tmagnitude = " << parameters.nearest_neighbor_strength << ";\n";
        of << "\tcutoff = " << parameters.nearest_neighbor_cutoff << ";\n";
        of << "\n};\n";

        of << "Dipole = {\n";
        of << "\tmagnitude = " << parameters.dipole_strength << ";\n";
        of << "\treciprocal mirrors = " << parameters.dipole_k_reciprocal << ";\n";
        of << "\treal mirrors = " << parameters.dipole_k_real << ";\n";
        of << "\n};";

        of << "RKKY = {\n";
        of << "\tmagnitude = " << parameters.rkky_strength << ";\n";
        of << "\tFermi wavevector = " << parameters.rkky_wavevector << ";\n";
        of << "\tcutoff = " << parameters.rkky_cutoff << ";\n";
        of << "};\n";

        of << "\n--- Geometry and size: ---\n\n";
        of << "Size = ";
        of << parameters.system_size[0] << " ";
        of << parameters.system_size[1] << " ";
        of << parameters.system_size[2] << ";\n";
    }
}