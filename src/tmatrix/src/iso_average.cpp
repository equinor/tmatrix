#include "TMatrix.hpp"

dTensor iso_average(dTensor const& arg)
{
    double lambda = (arg(0, 0) + arg(2, 2) + 5. * arg(0, 1) + 8. * arg(0, 2) -
                     2. * arg(3, 3)) /
                    15.;
    double mu = (7. * arg(0, 0) + 2. * arg(2, 2) - 5. * arg(0, 1) -
                 4. * arg(0, 2) + 6. * arg(3, 3)) /
                30.;

    double diag = lambda + 2. * mu;

    dTensor output = dTensor::Zero();
    output.topLeftCorner<3, 3>() << diag, lambda, lambda, lambda, diag, lambda,
        lambda, lambda, diag;
    output.bottomRightCorner<3, 3>().diagonal() =
        Eigen::Array3d::Constant(2. * mu);

    return output;
}
