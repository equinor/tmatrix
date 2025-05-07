#include "TMatrix.hpp"

dTensorVector calc_Kd(dTensor const& C0, dTensor const& C0_inv,
                      dTensorVector const& greenTensors)
{
    dTensorVector output;
    output.reserve(greenTensors.size());

    for (auto const& G : greenTensors) {
        output.emplace_back(dTensor::Zero());
        output.back().topLeftCorner<3, 3>() =
            (dSubTensor::Identity() +
             G.topLeftCorner<3, 3>() * C0.topLeftCorner<3, 3>())
                .inverse() *
            C0_inv.topLeftCorner<3, 3>();

        output.back().bottomRightCorner<3, 3>().diagonal() =
            C0_inv.bottomRightCorner<3, 3>().diagonal().array() /
            (1. + G.bottomRightCorner<3, 3>().diagonal().array() *
                      C0.bottomRightCorner<3, 3>().diagonal().array());
    }

    return output;
}
