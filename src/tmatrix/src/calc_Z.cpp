#include "TMatrix.hpp"
#include <cassert>

zTensorVector calc_Z(dTensor const& C0inv, dTensorVector const& td,
                     double omega, std::vector<double> const& gamma,
                     std::vector<double> const& v,
                     std::vector<double> const& tau)
{
    assert(td.size() == gamma.size());
    assert(gamma.size() == v.size());
    assert(v.size() == tau.size());

    zTensorVector output;
    output.reserve(td.size());

    zSubTensor reducer = zSubTensor::Zero();
    for (size_t i = 0; i < td.size(); i++)
        reducer += td[i].topLeftCorner<3, 3>() *
                   (v[i] / (1. + I * omega * gamma[i] * tau[i]));

    /* S0 * I2I2 * S0 */
    double val = std::pow(C0inv.topLeftCorner<1, 3>().sum(), 2);
    for (const auto & t : td) {
        output.emplace_back(zTensor::Zero());
        output.back().topLeftCorner<3, 3>() =
            t.topLeftCorner<3, 3>() * dSubTensor::Constant(val) * reducer;
    }
    // for (size_t i = 0; i < td.size(); i++) {
    //     output.emplace_back(zTensor::Zero());
    //     output.back().topLeftCorner<3, 3>() =
    //         td[i].topLeftCorner<3, 3>() * dSubTensor::Constant(val) * reducer;
    // }

    return output;
}
