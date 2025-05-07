#include "TMatrix.hpp"
#include <cassert>

dTensorVector iso_av_all(dTensorVector const& X, control1 ctrl)
{
    dTensorVector output;
    output.reserve(X.size());

    switch (ctrl) {
    case control1::ISOTROPIC:
            for (auto& tensor : X) output.emplace_back(iso_average(tensor));
            break;
    case control1::MIXED_CASE:
            assert(!(X.size() % 2));
            for (size_t i = 0; i < X.size() / 2; i++)
                output.emplace_back(iso_average(X[i]));
            for (size_t i = X.size() / 2; i < X.size(); i++)
                output.push_back(X[i]);
            break;
    case control1::ANISOTROPIC:
            return X;
    }

    return output;
}
