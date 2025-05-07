/* Regression test successful agains MATLAB */
#include "TMatrix.hpp"

dTensor Gtensor(dTensor const& C0, dTensor const& C0inv, double const alpha)
{
    /* The tensor C0 is block diagonal with (3,3) blocks, and diagonal in the
     * lower bottom block*/

    double mu = .5 * C0(3, 3);
    double kappa = C0(0, 0) - (4. / 3.) * mu;
    double poir = (1.5 * kappa - mu) / (3 * kappa + mu);

    dSubTensor S_r = dSubTensor();
    dTensor output = dTensor::Zero();

    if (alpha < 1) {
        double q =
            alpha / std::pow(1. - std::pow(alpha, 2), 1.5) *
            (std::acos(alpha) - alpha * std::sqrt(1. - std::pow(alpha, 2)));

        S_r(0, 0) = (3 / (8 * (1 - poir))) *
                        (std::pow(alpha, 2) / (std::pow(alpha, 2) - 1)) +
                    (1 / (4 * (1 - poir))) *
                        (1 - 2 * poir - 9 / (4 * (std::pow(alpha, 2) - 1))) * q;
        S_r(1, 0) = (1 / (4 * (1 - poir))) *
                    (std::pow(alpha, 2) / (2 * (std::pow(alpha, 2) - 1)) -
                     (1 - 2 * poir + 3 / (4 * (std::pow(alpha, 2) - 1))) * q);
        S_r(2, 0) = (1. / (2. * (1. - poir))) *
                    (-std::pow(alpha, 2) / (std::pow(alpha, 2) - 1.) +
                     0.5 *
                         (3. * std::pow(alpha, 2) / (std::pow(alpha, 2) - 1.) -
                          (1. - 2. * poir)) *
                         q);
        S_r(0, 1) = S_r(1, 0);
        S_r(1, 1) = S_r(0, 0);
        S_r(2, 1) = S_r(2, 0);
        S_r(0, 2) =
            (1. / (2. * (1. - poir))) *
            (2. * poir - 1. - 1. / (std::pow(alpha, 2) - 1.) +
             (1. - 2. * poir + 3. / (2. * (std::pow(alpha, 2) - 1.))) * q);
        S_r(1, 2) = S_r(0, 2);
        S_r(2, 2) = (1 / (2 * (1 - poir))) *
                    (1 - 2 * poir +
                     (3 * std::pow(alpha, 2) - 1) / (std::pow(alpha, 2) - 1) -
                     (1 - 2 * poir +
                      3 * std::pow(alpha, 2) / (std::pow(alpha, 2) - 1)) *
                         q);

        // Compute the lower right block
        Eigen::Array3d diag2;
        double S3 =
            (1 / (4 * (1 - poir))) *
            (1 - 2 * poir -
             (std::pow(alpha, 2) + 1) / (std::pow(alpha, 2) - 1) -
             0.5 *
                 (1 - 2 * poir -
                  (3 * (std::pow(alpha, 2) + 1)) / (std::pow(alpha, 2) - 1)) *
                 q);
        double S5 = (1 / (4 * (1 - poir))) *
                    (std::pow(alpha, 2) / (2 * (std::pow(alpha, 2) - 1)) +
                     (1 - 2 * poir - 3 / (4 * (std::pow(alpha, 2) - 1))) * q);
        diag2 << S3, S3, S5;

        Eigen::Array3d tmp = C0inv.bottomRightCorner<3, 3>().diagonal();
        output.bottomRightCorner<3, 3>().diagonal() = -2. * diag2 * tmp;
    } else if (1 == alpha) {
        double expr[3] = {(5. * poir - 1.) / (15. * (1. - poir)),
                          (2. * (4. - 5. * poir)) / (15. * (1. - poir)),
                          (4. - 5. * poir) / (15. * (1. - poir))};

        S_r(0, 0) = expr[0] + expr[1];
        S_r(1, 0) = expr[0];
        S_r(2, 0) = expr[0];
        S_r(0, 1) = expr[0];
        S_r(1, 1) = expr[0] + expr[1];
        S_r(2, 1) = expr[0];
        S_r(0, 2) = expr[0];
        S_r(1, 2) = expr[0];
        S_r(2, 2) = expr[0] + expr[1];

        // Compute the lower right block
        output.bottomRightCorner<3, 3>().diagonal() =
            -2. * expr[2] * C0inv.bottomRightCorner<3, 3>().diagonal();
    } else {
        return output;
    }

    output.topLeftCorner<3, 3>().noalias() =
        -(S_r.transpose() * (C0inv.topLeftCorner<3, 3>()));

    return output;
}
