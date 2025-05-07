#include "TMatrix.hpp"

dTensorVector calc_X(dTensor const& C0inv, dTensorVector const& td)
{
    dTensorVector output;
    output.reserve(td.size());
    dSubTensor const block = dSubTensor::Ones();

    for (auto const& obj : td) {
        output.emplace_back(dTensor::Zero());
        output.back().topLeftCorner<3, 3>() =
            obj.topLeftCorner<3, 3>() * C0inv.topLeftCorner<3, 3>() * block *
            C0inv.topLeftCorner<3, 3>() * obj.topLeftCorner<3, 3>();
    }

    return output;
}
