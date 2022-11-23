#include "include/couplings_kelvin.hpp"

namespace classmag::base
{
    MaterialData default_xas0(RareEarthIon ion)
    {
        MaterialData result;
        result.kondo_coupling_kelvin = 2321.2;
        result.length_scale_meters = 14.7e-10;
        result.atoms = 176.0;

        switch(ion)
        {
            case RareEarthIon::Gadolinium: 
                result.gJ = 7.0;
                result.G = 15.75;
                break;
            case RareEarthIon::Terbium:
                result.gJ = 9.0;
                result.G = 10.5;
                break;
            default:
                break;
        }

        return result;
    };

    double compute_dipole_magnitude(const MaterialData& parameters)
    {
        MagneticConstantsKelvin si;
        double result = pow(parameters.gJ,2.0) / 
            pow(parameters.length_scale_meters, 3.0);
        
        result *= si.mu_0 * pow(si.mu_b,2.0) /(4.0 * 3.1416 * si.k_B);

        return result;
    };

    double compute_rkky_magnitude(const MaterialData& params)
    {
        MagneticConstantsKelvin si;
        double rho = 
            params.atoms 
            * params.nu;
        
        double k_F = pow(3.0 * rho * pow(M_PI, 2.0),1.0/3.0);

        double result = 9.0 * M_PI/32.0;
        result *= pow(params.kondo_coupling_kelvin * si.k_B * params.nu/(k_F * si.hbar), 2.0);
        result *= pow(params.length_scale_meters, 2.0);
        result *= si.m_e * params.G/si.k_B;
        result *= pow(2.0 * k_F, -3.0);
        return result;
    };

    double compute_rkky_magnitude_kubler(const MaterialData& params)
    {
        MagneticConstantsKelvin si;
        double rho = 
            params.atoms 
            * params.nu *
            pow(params.length_scale_meters, -3.0);
        
        double k_F = pow(3.0 * pow(M_PI, 2.0) * rho,1.0/3.0);

        double result = pow(params.kondo_coupling_kelvin * si.k_B / si.hbar, 2.0);
        result *=  16.0 * si.m_e * pow(k_F, 4.0) * params.G;
        result *= 1.0/(si.k_B * pow(2.0 * M_PI, 3.0));
        return result;
    };

    CouplingsInKelvin compute_couplings(const MaterialData& parameters)
    {
        CouplingsInKelvin result;

        result.dipole_coupling = compute_dipole_magnitude(parameters);
        
        result.rkky_coupling = compute_rkky_magnitude(parameters);

        return result;
    };

    CouplingsInKelvin compute_xas0_couplings(RareEarthIon ion, double nu)
    {
        auto params = default_xas0(ion);
        params.nu = nu;

        auto result = compute_couplings(params);
        return result;
    };
}