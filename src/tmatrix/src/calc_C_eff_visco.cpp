#include "TMatrix.hpp"

#if 0
zTensor calc_C_eff_visco(
    double const velocity,                              // Vs
    double const permeability,                          // K_r
    double const viscocity,                             // eta_f
    std::vector<double> const& concentrationInclusion,  // v
    std::vector<double> const& gamma, std::vector<double> const& tau,
    std::vector<double> const& Kd_uuvv,
    double const materialBulkModulus,  // kappa
    double const fluidBulkModulus,     // kappa_f
    dTensor const& C0, dTensor const& C0inv, dTensor const& C1,
    dTensorVector const& tMatrices,  // td_bar
    dTensorVector const& Xtensor,    // X_bar
    dTensor const& corrMatrix,       // Gd
    double const frequency)
{
    double omega = 2 * M_PI * frequency;
    zdouble Theta =
        calc_Theta(concentrationInclusion, omega, gamma, tau, Kd_uuvv,
                   permeability / viscocity, omega / velocity,
                   materialBulkModulus, fluidBulkModulus);
    zTensorVector Z_bar =
        calc_Z(C0inv, tMatrices, omega, gamma, concentrationInclusion, tau);
    zTensorVector t_bar = calc_t(tMatrices, Theta, Xtensor, Z_bar, omega, gamma,
                                 tau, fluidBulkModulus);

    /*                        V                  V */
    return calc_C_eff(C0, C1, corrMatrix, t_bar, concentrationInclusion);
}

#else

// Optimized version where the calc_Theta, calc_Z and calc_t functions are
// inlined, and all overhead and superfluous allocations are removed. Tests 
// indicate around 30% performance gain 
zTensor calc_C_eff_visco(double velocity, double permeability, double viscocity,
                         const std::vector<double>& concentrationInclusion,
                         const std::vector<double>& gamma,
                         const std::vector<double>& tau,
                         const std::vector<double>& Kd_uuvv,
                         double materialBulkModulus, double fluidBulkModulus,
                         const dTensor& C0, const dTensor& C0inv,
                         dTensor const& C1, const dTensorVector& tMatrices,
                         const dTensorVector& Xtensor,
                         const dTensor& corrMatrix, double frequency)
{
    const double omega = 2. * M_PI * frequency;
    const zdouble Theta =
        calc_Theta(concentrationInclusion, omega, gamma, tau, Kd_uuvv,
                   permeability / viscocity, omega / velocity,
                   materialBulkModulus, fluidBulkModulus);

    zSubTensor reducer = zSubTensor::Zero();
    for (size_t i = 0; i < tMatrices.size(); i++)
        reducer +=
            tMatrices[i].topLeftCorner<3, 3>() *
            (concentrationInclusion[i] / (1. + I * omega * gamma[i] * tau[i]));

    const double val = std::pow(C0inv.topLeftCorner<1, 3>().sum(), 2);

    zTensor C1cpy = C1.cast<zdouble>();

    for (size_t i = 0; i < tMatrices.size(); i++) {
        zSubTensor Z_upper = tMatrices[i].topLeftCorner<3, 3>() *
                       dSubTensor::Constant(val) * reducer;

        zSubTensor t = tMatrices[i].topLeftCorner<3, 3>().cast<zdouble>() +
            (Theta * Z_upper +
             I * omega * tau[i] * fluidBulkModulus *
                 Xtensor[i].topLeftCorner<3, 3>().cast<zdouble>()) /
                (1. + I * omega * gamma[i] * tau[i]);

        for (Eigen::Index j = 3; j < 6; j++) {
            C1cpy(j, j) +=
                concentrationInclusion[i] *
                (tMatrices[i](j, j) +
                 I * omega * tau[i] * fluidBulkModulus * Xtensor[i](j, j) /
                     ((1. + I * omega * gamma[i] * tau[i])));
        }
        C1cpy.topLeftCorner<3, 3>() += concentrationInclusion[i] * t;
    }

    zSubTensor tmp1 = (zSubTensor::Identity() +
                 corrMatrix.topLeftCorner<3, 3>().cast<zdouble>() *
                     C1cpy.topLeftCorner<3, 3>())
                    .inverse();

    zTensor tmp2 = zTensor::Zero();
    tmp2.topLeftCorner<3, 3>().noalias() = C1cpy.topLeftCorner<3, 3>() * tmp1;
    tmp2.bottomRightCorner<3, 3>().diagonal() =
        lowDiag(C1cpy).array() /
        (1. +
         lowDiag(corrMatrix).cast<zdouble>().array() * lowDiag(C1cpy).array());

    return C0.cast<zdouble>() + tmp2;
}
#endif
