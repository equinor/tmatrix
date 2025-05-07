#pragma once

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

const double PERMEABILITY_SU_UNIT_CONV = 0.986923e-15;
const double FLUID_SU_UNIT_CONV = 1e-2;

typedef std::complex<double> zdouble;
typedef Eigen::Matrix<double, 6, 6> dTensor;
typedef Eigen::Matrix<double, 3, 3> dSubTensor;
typedef Eigen::Matrix<zdouble, 6, 6> zTensor;
typedef Eigen::Matrix<zdouble, 3, 3> zSubTensor;
typedef std::vector<dTensor, Eigen::aligned_allocator<dTensor>> dTensorVector;
typedef std::vector<zTensor, Eigen::aligned_allocator<zTensor>> zTensorVector;

const zdouble I(0.0, 1.0);

#include "TMatrixAPI.hpp"

static inline auto lowDiag(dTensor const& in_)
    -> decltype(in_.bottomRightCorner<3, 3>().diagonal())
{
    return in_.bottomRightCorner<3, 3>().diagonal();
}

static inline auto lowDiag(zTensor const& in_)
    -> decltype(in_.bottomRightCorner<3, 3>().diagonal())
{
    return in_.bottomRightCorner<3, 3>().diagonal();
}

dTensor Gtensor(dTensor const& C0, dTensor const& C0inv, double const alpha);

zTensorVector calc_t(dTensorVector const& td, zdouble const Theta,
                     dTensorVector const& X, zTensorVector const& Z,
                     double const omega, std::vector<double> const& gamma,
                     std::vector<double> const& tau,
                     double const fluidBulkModulus);

dTensorVector calc_td(dTensor const& C0, dTensorVector const& greenTensors,
                      dTensorVector const& Kd);

dTensor calc_isolated(dTensor const& C0, double const fluidBulkModulus,
                      dTensorVector const& greenTensors,
                      std::vector<double> const& concentrations,
                      control1 ctrl);

dTensorVector iso_av_all(dTensorVector const& X, control1 ctrl);

dTensorVector calc_X(dTensor const& C0inv, dTensorVector const& td);

dTensorVector calc_Kd(dTensor const& C0, dTensor const& C0_inv,
                      dTensorVector const& greenTensors);

zTensorVector calc_Z(dTensor const& C0inv, dTensorVector const& td,
                     double omega, std::vector<double> const& gamma,
                     std::vector<double> const& v,
                     std::vector<double> const& tau);

std::vector<double> calc_Kd_uuvv(dTensorVector const& Kd);

zTensor calc_C_eff(dTensor const& C0, dTensor const& C1, dTensor const& Gd,
                   zTensorVector const& t, std::vector<double> v);

dTensor iso_average(dTensor const& arg);

Result velocity_vti_angles(zTensor const& Ceff, double const rho,
                           double const rad_angle);

zdouble calc_Theta(std::vector<double> cavityConcentration, double omega,
                   std::vector<double> gamma, std::vector<double> tau,
                   std::vector<double> Kd_uuvv, double permViscRatio,
                   double wavenumber, double materialBulkModulus,
                   double fluidBulkModulus);

zTensor calc_C_eff_visco(
    double const velocity, double const permeability, double const viscocity,
    std::vector<double> const& concentrationInclusion,
    std::vector<double> const& gamma, std::vector<double> const& tau,
    std::vector<double> const& Kd_uuvv, double const materialBulkModulus,
    double const fluidBulkModulus, dTensor const& C0, dTensor const& C0inv,
    dTensor const& C1, dTensorVector const& tMatrices,
    dTensorVector const& Xtensor, dTensor const& corrMatrix,
    double const frequency);
