#include "TMatrix.hpp"
#include <cassert>

zTensorVector calc_t(dTensorVector const& td, zdouble const Theta,
                     dTensorVector const& X, zTensorVector const& Z,
                     double const omega, std::vector<double> const& gamma,
                     std::vector<double> const& tau,
                     double const fluidBulkModulus)
{
    assert(td.size() == gamma.size());
    assert(gamma.size() == tau.size());

    zTensorVector output;
    output.reserve(td.size());

    zTensor t = zTensor::Zero();
    for (size_t i = 0; i < td.size(); i++) {
        t.topLeftCorner<3, 3>() =
            td[i].topLeftCorner<3, 3>().cast<zdouble>() +
            (Theta * Z[i].topLeftCorner<3, 3>() +
             I * omega * tau[i] * fluidBulkModulus *
                 X[i].topLeftCorner<3, 3>().cast<zdouble>()) /
                (1. + I * omega * gamma[i] * tau[i]);

        for (Eigen::Index j = 3; j < 6; j++)
            t(j, j) = td[i](j, j) + Theta * Z[i](j, j) +
                      I * omega * tau[i] * fluidBulkModulus * X[i](j, j) /
                          ((1. + I * omega * gamma[i] * tau[i]));

        /*
                t.bottomRightCorner<3, 3>().diagonal() =
                    lowDiag(td[i]).array().cast<zdouble>() +
                    (Theta * lowDiag(Z[i]).array() +
                     I * omega * tau[i] * fluidBulkModulus *
                         lowDiag(X[i]).array().cast<zdouble>()) /
                        (1. + I * omega * gamma[i] * tau[i]);

        */

        output.push_back(t);
    }

    return output;
}
