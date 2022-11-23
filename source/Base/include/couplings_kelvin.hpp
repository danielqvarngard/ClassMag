#ifndef CLASSMAG_BASE_SIVALUES_HPP
#define CLASSMAG_BASE_SIVALUES_HPP

#include "math.h"
#include "cmath"
namespace classmag::base
{
    struct MagneticConstantsKelvin
    {
        const double hbar = 1.0546e-34;
        const double m_e = 9.109e-31;
        const double mu_b = 9.274e-24;
        const double mu_0 = 4.0 * 3.1416 * 1.0e-7;
        const double k_B = 1.387e-23;
    };

    enum class RareEarthIon
    {
        Gadolinium,
        Terbium
    };

    struct CouplingsInKelvin
    {
        double dipole_coupling;
        double rkky_coupling;
    };

    struct MaterialData
    {
        double length_scale_meters;
        double kondo_coupling_kelvin;
        double G;
        double gJ;
        double nu; // electrons per atom for Fermi wavevector computation
        double atoms; // total atoms per unit cell
    };

    MaterialData default_xas0(RareEarthIon ion);

    CouplingsInKelvin compute_couplings(const MaterialData& parameters);

    double compute_dipole_magnitude(const MaterialData& parameters);

    double compute_rkky_magnitude(const MaterialData& parameters);

    double compute_rkky_magnitude_kubler(const MaterialData& parameters);

    CouplingsInKelvin compute_xas0_couplings(RareEarthIon ion, double nu);
}


#endif