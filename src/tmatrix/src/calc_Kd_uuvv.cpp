#include "TMatrix.hpp"

std::vector<double> calc_Kd_uuvv(dTensorVector const& Kd)
{
    std::vector<double> output;
    output.reserve(Kd.size());

    for (auto const& tensor : Kd)
        output.push_back(tensor.topLeftCorner<3, 3>().sum());

    return output;
}
