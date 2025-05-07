#include <iostream>
#include <vector>
#include <array>
#include <cstddef>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "TMatrixAPI.hpp"

namespace py = pybind11;

int tmatrix_porosity(
    py::array_t<double> out_np,
    size_t dim,
    py::array_t<double> mineral_property_np,
    py::array_t<double> fluid_property_np,
    py::array_t<double> phi_vector_np,
	int in_scenario,
    double frequency,
    double angle_of_sym_plane,
    double per_inc_con,
    double per_inc_ani
) {
    auto out = out_np.mutable_data();
    auto mineral_property = mineral_property_np.data();
    auto fluid_property = fluid_property_np.data();
    auto phi_vector = phi_vector_np.data();

	scenario scen;
	switch (in_scenario) {
	case 1:
		scen = scenario::DUAL_POR_MOSTLY_ROUNDED;
		break;
	case 2:
		scen = scenario::DUAL_POR_LITTLE_ROUNDED;
		break;
	case 3:
		scen = scenario::MIXED_PORES;
		break;
	case 4:
		scen = scenario::FLAT_PORES_AND_CRACKS;
		break;
	default:
		scen = scenario::DUAL_POR_MOSTLY_ROUNDED;
	}

	TMatrix_Porosity t_matrix(scen, per_inc_con, per_inc_ani);

	std::array<double, 3> mineral;
	std::array<double, 4> fluid;

	for (size_t n = 0; n < dim; n++)
	{
		std::memcpy(mineral.data(), &mineral_property[n * 3], 3*sizeof(double));
		std::memcpy(fluid.data(), &fluid_property[n * 4], 4*sizeof(double));

		Result res = t_matrix.evaluate(mineral, fluid, phi_vector[n], frequency, angle_of_sym_plane);

		double* out_n = &out[n * 4];
		out_n[0] = res.pVelocity;
		out_n[1] = res.horzShearVelocity;
		out_n[2] = res.vertShearVelocity;
		out_n[3] = res.rhoEffective;
	}
	return 0;
}

void tmatrix_porosity_noscenario(
    py::array_t<double> out_np,
    size_t out_N,
    py::array_t<double> mineral_property_np,
    py::array_t<double> fluid_property_np,
    py::array_t<double> phi_vector_np,
    py::array_t<double> alpha_np,
    py::array_t<double> v_np,
    py::array_t<int> alpha_size_np,
    size_t alpha_N,
    double frequency,
    double angle,
    py::array_t<double> inc_con_np,
    py::array_t<double> inc_ani_np,
    size_t inc_con_N
) {
    auto out = out_np.mutable_data();

    auto mineral_property = mineral_property_np.data();
    auto fluid_property = fluid_property_np.data();

    auto phi_vector = phi_vector_np.data();
    auto alpha = alpha_np.data();
    auto v = v_np.data();
    auto alpha_size = alpha_size_np.data();
    auto inc_con = inc_con_np.data();
    auto inc_ani = inc_ani_np.data();

	// First sample
	TMatrix_Porosity t_matrix(inc_con[0], inc_ani[0], alpha_size[0]);

	// Read in first alpha and V sample
	std::vector<double> alpha_v(alpha, alpha + alpha_size[0]);
	std::vector<double> v_v(v, v + alpha_size[0]);
	t_matrix.setAlpha(alpha_v);
	t_matrix.setV(v_v);

	std::array<double, 3> mineral;
	std::array<double, 4> fluid;
	if ((alpha_N == 1) && (inc_con_N == 1))
	{
		for (size_t n = 0; n < out_N; n++)
		{
			std::memcpy(mineral.data(), &mineral_property[n * 3], 3 * sizeof(double));
			std::memcpy(fluid.data(), &fluid_property[n * 4], 4 * sizeof(double));

			Result res = t_matrix.evaluate(mineral, fluid, phi_vector[n], frequency, angle);

			double* out_n = &out[n * 4];
			out_n[0] = res.pVelocity;
			out_n[1] = res.horzShearVelocity;
			out_n[2] = res.vertShearVelocity;
			out_n[3] = res.rhoEffective;
		}
	}
	else if ((alpha_N == out_N) && (inc_con_N == 1))
	{
		size_t alpha_id_start = 0;
		for (size_t n = 0; n < out_N; n++)
		{
			std::memcpy(mineral.data(), &mineral_property[n * 3], 3 * sizeof(double));
			std::memcpy(fluid.data(), &fluid_property[n * 4], 4 * sizeof(double));

			alpha_v = std::vector<double>(alpha + alpha_id_start, alpha + alpha_id_start + alpha_size[n]);
			v_v = std::vector<double>(v + alpha_id_start, v + alpha_id_start + alpha_size[n]);
			alpha_id_start += alpha_size[n];

			Result res = t_matrix.evaluate(alpha_v, v_v, mineral, fluid, phi_vector[n], frequency, angle);

			double* out_n = &out[n * 4];
			out_n[0] = res.pVelocity;
			out_n[1] = res.horzShearVelocity;
			out_n[2] = res.vertShearVelocity;
			out_n[3] = res.rhoEffective;
		}
	}
	else if ((alpha_N == out_N) && (inc_con_N == out_N))
	{
		size_t alpha_id_start = 0;
		for (size_t n = 0; n < out_N; n++)
		{
			std::memcpy(mineral.data(), &mineral_property[n * 3], 3 * sizeof(double));
			std::memcpy(fluid.data(), &fluid_property[n * 4], 4 * sizeof(double));

			alpha_v = std::vector<double>(alpha + alpha_id_start, alpha + alpha_id_start + alpha_size[n]);
			v_v = std::vector<double>(v + alpha_id_start, v + alpha_id_start + alpha_size[n]);
			alpha_id_start += alpha_size[n];


			Result res = t_matrix.evaluate(inc_con[n], inc_ani[n], alpha_v, v_v,
				mineral, fluid, phi_vector[n], frequency, angle);

			double* out_n = &out[n * 4];

			out_n[0] = res.pVelocity;
			out_n[1] = res.horzShearVelocity;
			out_n[2] = res.vertShearVelocity;
			out_n[3] = res.rhoEffective;
		}
	}
}


PYBIND11_MODULE(_tmatrix, m) {
    m.def("tmatrix_porosity", &tmatrix_porosity, "Compute TMatrix Porosity",
          py::arg("out_np"),
          py::arg("dim"),
          py::arg("mineral_property_np"),
          py::arg("fluid_property_np"),
          py::arg("phi_vector_np"),
          py::arg("in_scenario"),
          py::arg("frequency"),
          py::arg("angle_of_sym_plane"),
          py::arg("per_inc_con"),
          py::arg("per_inc_any"));

    m.def("tmatrix_porosity_noscenario", &tmatrix_porosity_noscenario,
            "Compute TMatrix Porosity, No scenario",
            py::arg("out_np"),
            py::arg("out_N"),
            py::arg("mineral_property_np"),
            py::arg("fluid_property_np"),
            py::arg("phi_vector_np"),
            py::arg("alpha_np"),
            py::arg("v_np"),
            py::arg("alpha_size_np"),
            py::arg("alpha_N"),
            py::arg("frequency"),
            py::arg("angle"),
            py::arg("inc_con_np"),
            py::arg("inc_ani_np"),
            py::arg("inc_con_N"));
}