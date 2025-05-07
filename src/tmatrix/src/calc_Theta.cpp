#include "TMatrix.hpp"

/* Note that in the original documentation it is stated that the function
 * returns a tensor, yet its code returs a scalar, and the original article (2nd
 * article in Agerborg's thesis, p.6 eq.(11)) also seems to indicate a scalar on
 * the lhs.
 */
zdouble calc_Theta(std::vector<double> cavityConcentration,  // v
                   double omega, std::vector<double> gamma,
                   std::vector<double> tau, std::vector<double> Kd_uuvv,
                   double permViscRatio,        // Dr
                   double wavenumber,           // k
                   double materialBulkModulus,  // kappa
                   double fluidBulkModulus)     // kappa_f
{
    assert(cavityConcentration.size() == gamma.size());
    assert(cavityConcentration.size() == tau.size());
    assert(cavityConcentration.size() == Kd_uuvv.size());

    zdouble sigma_a = 0;
    zdouble sigma_b = 0;
    for (size_t i = 0; i < cavityConcentration.size(); i++) {
        zdouble tmp =
            cavityConcentration[i] / (1. + I * omega * gamma[i] * tau[i]);
        sigma_a += tmp;
        sigma_b += tmp * Kd_uuvv[i];
    }

    return fluidBulkModulus /
           ((1. - fluidBulkModulus / materialBulkModulus) * sigma_a +
            fluidBulkModulus * sigma_b -
            (I * std::pow(wavenumber, 2) / omega) * permViscRatio *
                fluidBulkModulus);
}
