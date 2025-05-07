#include "TMatrix.hpp"
#include <cassert>

dTensorVector calc_td(dTensor const& C0, dTensorVector const& greenTensors,
                      dTensorVector const& Kd)
{
    assert(greenTensors.size() == Kd.size());

    dTensorVector output;
    output.reserve(Kd.size());

    size_t i = 0;
    Eigen::Array3d const lowerC0_diag = C0.bottomRightCorner<3, 3>().diagonal();

    for (auto const& Ktensor : Kd) {
        dTensor const& G = greenTensors[i++];
        output.emplace_back(dTensor::Zero());

        output.back().topLeftCorner<3, 3>().noalias() =
            G.topLeftCorner<3, 3>().inverse() *
            (Ktensor.topLeftCorner<3, 3>() * C0.topLeftCorner<3, 3>() -
             dSubTensor::Identity());

        Eigen::Array3d tmp = lowDiag(Ktensor).array() * lowerC0_diag;
        output.back().bottomRightCorner<3, 3>().diagonal() =
            (tmp - 1.) / lowDiag(G).array();
    }

    return output;
}
