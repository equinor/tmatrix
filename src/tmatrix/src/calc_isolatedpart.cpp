#include "TMatrix.hpp"
#include <cassert>

void correctionTensorLoop(dTensor const& G, dTensor const& C0,
                          double const /*fluidBulkModulus*/, dTensor& output,
                          control1 ctrl,
                          std::vector<double> const& concentrations, size_t i);

dTensor calc_isolated(dTensor const& C0, double const fluidBulkModulus,
                      dTensorVector const& greenTensors,
                      std::vector<double> const& concentrations,
                      control1 ctrl)
{
    dTensor output = dTensor::Zero();
    size_t i = 0;

    auto loopfun = [&](dTensor const& G_, control1 ctrl_) -> void {
        return correctionTensorLoop(G_, C0, fluidBulkModulus, output, ctrl_,
                                    concentrations, i++);
    };

    if (ctrl != control1::MIXED_CASE) {
        for (auto const& G : greenTensors) loopfun(G, ctrl);
    } else {  // MIXED_CASE ctrl
        assert(!(greenTensors.size() % 2));

        for (size_t j = 0; j < greenTensors.size() / 2; j++)
            loopfun(greenTensors[j], control1::ISOTROPIC);

        for (size_t j = greenTensors.size() / 2; j < greenTensors.size(); j++)
            loopfun(greenTensors[j], control1::ANISOTROPIC);
    }

    return output;
}

void correctionTensorLoop(dTensor const& G, dTensor const& C0,
                          double const fluidBulkModulus, dTensor& output,
                          control1 ctrl,
                          std::vector<double> const& concentrations, size_t i)
{
    dTensor tmp = dTensor::Zero();

    // For some reason Eigen and MATLAB start diverging on the 10th digit in
    // this expression...
    tmp.topLeftCorner<3, 3>() =
        (dSubTensor::Constant(fluidBulkModulus) - C0.topLeftCorner<3, 3>()) *
        (dSubTensor::Identity() -
         G.topLeftCorner<3, 3>() * (dSubTensor::Constant(fluidBulkModulus) -
                                    C0.topLeftCorner<3, 3>()))
            .inverse();

    tmp.bottomRightCorner<3, 3>().diagonal() =
        -lowDiag(C0).array() / (1. + lowDiag(G).array() * lowDiag(C0).array());

    if (ctrl != control1::ANISOTROPIC) tmp = iso_average(tmp);

    output.topLeftCorner<3, 3>() +=
        concentrations[i] * tmp.topLeftCorner<3, 3>();
    output.bottomRightCorner<3, 3>().diagonal() +=
        concentrations[i++] * lowDiag(tmp);
}
