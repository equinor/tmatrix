#include "TMatrix.hpp"
#include <algorithm>
#include <iterator>
#include <vector>

namespace {

template <typename T>
std::vector<T> concat(const std::vector<T>& A, const std::vector<T>& B)
{
    std::vector<T> AB;
    AB.reserve(A.size() + B.size());          // to avoid reallocation
    AB.insert(AB.end(), A.begin(), A.end());  // insert A before end
    AB.insert(AB.end(), B.begin(), B.end());  // insert B before end
    return AB;
}

}  // namespace

using dlist = std::initializer_list<double>;
static const std::initializer_list<double> static_alphas[] = {
    dlist{.9, .1},             // DUAL_POR_MOSTLY_ROUNDED
    dlist{.58, .027},          // DUAL_POR_LITTLE_ROUNDED
    dlist{.9, .1, 1e-2},       // MIXED_PORES
    dlist{.9, .1, 1e-2, 1e-3}  // FLAT_PORES_AND_CRACK
};

static const std::initializer_list<double> static_vs[] = {
    dlist{.9, .1},               // DUAL_POR_MOSTLY_ROUNDED
    dlist{.85, .15},             // DUAL_POR_LITTLE_ROUNDED
    dlist{.8, .19, 1e-2},        // MIXED_PORES
    dlist{.689, .3, 1e-2, 1e-3}  // FLAT_PORES_AND_CRACK
};

TMatrix_Porosity::TMatrix_Porosity(double per_inc_con, double per_inc_ani,
                                   size_t len)
    : tau(std::vector<double>(len, 1e-7)),
      per_inc_con(per_inc_con),
      per_inc_ani(per_inc_ani)
{
    updateControl();
    if (per_inc_ani != 0) tau = concat(tau, tau);
}

void TMatrix_Porosity::updateControl()
{
    if (per_inc_ani != 0. && per_inc_ani != 1.)
        ctrl1 = control1::MIXED_CASE;
    else if (per_inc_ani == 1.)
        ctrl1 = control1::ANISOTROPIC;
    else
        ctrl1 = control1::ISOTROPIC;

    if (per_inc_con != 0. && per_inc_con != 1.)
        ctrl2 = control2::BOTH;
    else if (per_inc_con == 1)
        ctrl2 = control2::CONNECTED;
    else
        ctrl2 = control2::ISOLATED;
}

// With constant alpha, the alpha goes out of scope before evaluate is called.
// Therefore we need a copy. As it is a small vector we take the overhead for
// the other cases
void TMatrix_Porosity::setAlpha(std::vector<double> new_alpha)
{
    alpha = std::move(new_alpha);

    size_t alen = alpha.size();
    if (per_inc_ani != 0. && per_inc_ani != 1.) alen *= 2;

    tau = std::vector<double>(alen, 1e-7);

    if (per_inc_ani != 0 && per_inc_ani != 1) {
        // append the vector to itself
        alpha = concat(alpha, alpha);
    }
}

void TMatrix_Porosity::setV(std::vector<double> new_v)
{
    v = std::move(new_v);
    if (per_inc_ani != 0 && per_inc_ani != 0) {
        v = concat(v, v);
        for (size_t i = 0; i < v.size() / 2; i++) v[i] *= (1. - per_inc_ani);
        for (size_t i = v.size() / 2; i < v.size(); i++) v[i] *= per_inc_ani;
    }
}

TMatrix_Porosity::TMatrix_Porosity(scenario scen, double per_inc_con,
                                   double per_inc_ani)

    : TMatrix_Porosity(
          per_inc_con, per_inc_ani,
          static_alphas[static_cast<std::underlying_type<scenario>::type>(scen)]
              .size())
{
    auto index = static_cast<std::underlying_type<scenario>::type>(scen);
    setAlpha(static_alphas[index]);
    setV(static_vs[index]);
}

Result TMatrix_Porosity::evaluate(std::array<double, 3> const& mineralProperty,
                                  std::array<double, 4> const& fluidProperty,
                                  double const porosity, double const frequency,
                                  double const angle) const
{
    if (frequency <= 0.)
        throw std::runtime_error("T-Matrix frequency must be positive");

    std::vector<double> in_v_connected;
    std::vector<double> in_v_isolated;

    for (auto& obj : v) in_v_connected.push_back(obj * per_inc_con * porosity);
    for (auto& obj : v)
        in_v_isolated.push_back(obj * (1. - per_inc_con) * porosity);

    for (size_t i = 0; i < in_v_connected.size(); i++)
        if (in_v_connected[i] * porosity > alpha[i])
            in_v_connected[i] = alpha[i] / 2.;

    for (size_t i = 0; i < in_v_isolated.size(); i++)
        if (in_v_isolated[i] * porosity > alpha[i])
            in_v_isolated[i] = alpha[i] / 2.;

    /* Move to class later */
    double kappa_f = fluidProperty[0];
    double rho_f = fluidProperty[1];
    double K_r = fluidProperty[2] * PERMEABILITY_SU_UNIT_CONV;
    double eta_f = fluidProperty[3] * FLUID_SU_UNIT_CONV;

    double rho_s = mineralProperty[2];
    double Vs_s = std::sqrt(mineralProperty[1] / rho_s);
    double Vp_s = std::sqrt(
        (mineralProperty[0] + (4. / 3.) * mineralProperty[1]) / rho_s);

    /*verified */
    dTensor C0 = dTensor::Zero();

    double expr[2] = {rho_s * std::pow(Vp_s, 2),
                      rho_s * (std::pow(Vp_s, 2) - 2. * std::pow(Vs_s, 2))};

    C0.topLeftCorner<3, 3>() = dSubTensor::Constant(expr[1]);
    C0.topLeftCorner<3, 3>().diagonal() = Eigen::Array3d::Constant(expr[0]);
    C0.bottomRightCorner<3, 3>().diagonal() =
        Eigen::Array3d::Constant(2. * rho_s * std::pow(Vs_s, 2));

    dTensor C0_inv = C0.inverse();

    /* verified */
    double kappa = C0(0, 0) - (4. / 3.) * rho_s * std::pow(Vs_s, 2);

    dTensor Gd = Gtensor(C0, C0_inv, 1.);
    dTensor C1;
    zTensor C_eff;

    if (ctrl2 != control2::CONNECTED) {
        dTensorVector G;
        for (auto& a : alpha) /*isolated_alpha */
            G.emplace_back(Gtensor(C0, C0_inv, a));

        C1 = calc_isolated(C0, kappa_f, G, in_v_isolated, ctrl1);
    } else
        C1 = dTensor::Zero();

    if (ctrl2 != control2::ISOLATED) {
        dTensorVector G;
        for (auto& a : alpha) /* connected_alpha */
            G.emplace_back(Gtensor(C0, C0_inv, a));

        dTensorVector Kd = calc_Kd(C0, C0_inv, G);

        dTensorVector td = calc_td(C0, G, Kd);
        dTensorVector td_bar = iso_av_all(td, ctrl1);
        dTensorVector X = calc_X(C0_inv, td);
        dTensorVector X_bar = iso_av_all(X, ctrl1);

        std::vector<double> Kd_uuvv = calc_Kd_uuvv(Kd);

        std::vector<double> gamma;
        gamma.reserve(Kd_uuvv.size());
        for (auto k : Kd_uuvv)
            gamma.push_back(1. - kappa_f / kappa + kappa_f * k);
        // for (size_t i = 0; i < Kd_uuvv.size(); i++)
        //     gamma.push_back(1. - kappa_f / kappa + kappa_f * Kd_uuvv[i]);

        C_eff = calc_C_eff_visco(Vs_s, K_r, eta_f, in_v_connected, gamma, tau,
                                 Kd_uuvv, kappa, kappa_f, C0, C0_inv, C1,
                                 td_bar, X_bar, Gd, frequency);
    } else
        /* Must be optimized for the matrix structure */
        C_eff = (C0 + C1 * (dTensor::Identity() + Gd * C1).inverse())
                    .cast<zdouble>();

    double phi = 0;
    if (ctrl2 != control2::ISOLATED)
        for (auto& obj : in_v_connected) phi += obj;

    if (ctrl2 != control2::CONNECTED)
        for (auto& obj : in_v_isolated) phi += obj;

    double rhoEffective = phi * rho_f + (1. - phi) * rho_s;

    Result res = velocity_vti_angles(
        C_eff, rhoEffective, ctrl1 == control1::ISOTROPIC ? 0 : angle * M_PI / 180.);
    res.frequency = frequency;
    res.rhoEffective = rhoEffective;

    return res;
}
