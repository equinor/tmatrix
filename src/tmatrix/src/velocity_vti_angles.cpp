#include "TMatrix.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <complex>

/* Returns the P-velocity and S-velocities */
Result velocity_vti_angles(zTensor const& Ceff, double const rho,
                           double const rad_angle)
{
    auto CC = Ceff.diagonal();

    /* {Vp90, Vp0, Vs0, Vs90} */
    Eigen::Array4d diag;
    diag << std::real(1. / std::sqrt(CC(0) / rho)),
        std::real(1. / std::sqrt(CC(2) / rho)),
        std::real(1. / std::sqrt(CC(3) / (2. * rho))),
        std::real(1. / std::sqrt(CC(5) / (2. * rho)));
    diag = 1. / diag;

    Eigen::Array4d Ccoefs = rho * diag.square();
    const double& C11 = Ccoefs[0];
    const double& C33 = Ccoefs[1];
    const double& C44 = Ccoefs[2];
    const double& C66 = Ccoefs[3];

    zdouble m45 =
        std::pow(.5 * (CC(0) - (.5 * CC(3))) - .5 * (CC(2) - .5 * CC(3)), 2) +
        std::pow(CC(0) - .5 * CC(3), 2);

    double Vp45 = 1. / std::real(1. / std::sqrt((.5 * (CC(0) + CC(2) + CC(3)) +
                                                 std::sqrt(m45)) /
                                                (2. * rho)));

    double C13 =
        -C44 + std::sqrt(4. * std::pow(rho, 2) * std::pow(Vp45, 4) -
                         2. * rho * std::pow(Vp45, 2) * (C11 + C33 + 2. * C44) +
                         (C11 + C44) * (C33 + C44));
    /*     m_real =         ((diag[0] -diag[2]) * (        sin(
     * rad_angle)^2)  - (diag[1] - diag[2]) * (        cos(     rad_angle)^2))^2
     * + (       (C13 + diag[2])^2) * (        sin(     2  * rad_angle)^2 );*/

    double sin_rad_sq = std::pow(std::sin(rad_angle), 2);
    double cos_rad_sq = std::pow(std::cos(rad_angle), 2);

    double m_real =
        std::pow((C11 - C44) * sin_rad_sq - (C33 - C44) * cos_rad_sq, 2) +
        std::pow(C13 + C44, 2) * std::pow(std::sin(2. * rad_angle), 2);

    Result res;
    res.pVelocity = std::sqrt(C11 * sin_rad_sq + C33 * cos_rad_sq + C44 +
                              std::sqrt(m_real)) *
                    std::sqrt(1. / (2. * rho));

    // TODO: horizontal and vertical shear are swapped with respect to original
    // code. But fixing breaks all the regression tests :-(

    res.horzShearVelocity = std::sqrt(C11 * sin_rad_sq + C33 * cos_rad_sq +
                                      C44 - std::sqrt(m_real)) *
                            std::sqrt(1. / (2. * rho));

    res.vertShearVelocity =
        std::sqrt((C66 * sin_rad_sq + C44 * cos_rad_sq) / rho);

    res.gamma = (C66 - C44) / (2. * C44);
    res.delta = (std::pow(C13 + C44, 2) - std::pow(C33 - C44, 2)) /
                (2. * C33 * (C33 - C44));
    res.epsilon = (C11 - C33) / (2. * C33);

    return res;
}
