#include "TMatrix.hpp"
#include <cassert>

zTensor calc_C_eff(dTensor const& C0, dTensor const& C1, dTensor const& Gd,
                   zTensorVector const& t, std::vector<double> v)
{
    zTensor C1cpy = C1.cast<zdouble>();
    assert(t.size() == v.size());

    for (size_t i = 0; i < t.size(); i++) C1cpy += v[i] * t[i];

    zSubTensor tmp1 =
        (zSubTensor::Identity() +
         Gd.topLeftCorner<3, 3>().cast<zdouble>() * C1cpy.topLeftCorner<3, 3>())
            .inverse();

    zTensor tmp2 = zTensor::Zero();
    tmp2.topLeftCorner<3, 3>().noalias() = C1cpy.topLeftCorner<3, 3>() * tmp1;
    tmp2.bottomRightCorner<3, 3>().diagonal() =
        lowDiag(C1cpy).array() /
        (1. + lowDiag(Gd).cast<zdouble>().array() * lowDiag(C1cpy).array());

    return C0.cast<zdouble>() + tmp2;
}
