#ifndef CLASSMAG_BASE_SIVALUES_HPP
#define CLASSMAG_BASE_SIVALUES_HPP

namespace classmag::base
{
    struct MagneticConstantsKelvin
    {
        const double hbar = 1.0546e-34;
        const double m_e = 9.109e-31;
        const double mu_b = 9.274e-24;
        const double mu_0 = 4.0 * 3.1416 * 10^-7;
        const double k_B = 1.387e-23;
    };

    enum class RareEarthIon
    {
        Gadolinium,
        Terbium
    };

    struct CouplingsInKelvin
    {
        double dipole_coupling,
        double rkky_coupling
    };

    struct MaterialData
    {
        double length_scale_meters,
        double G,
        double gJ
    };

    MaterialData default_xas0(RareEarthIon ion);
}


#endif