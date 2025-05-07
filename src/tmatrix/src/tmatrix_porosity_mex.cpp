#include "mex.h"
#include "TMatrixAPI.hpp"
#include <array>
#include <functional>
#include <string>
#include <tuple>
#include <cassert>

#define N_OUT_FLD 5

const size_t SCALAR = 1;

const char* fieldnames[N_OUT_FLD] = {"Vp", "Vsv", "Vsh", "frequency",
                                     "rho_eff"};

using tmatrix_eval = std::function<Result(TMatrix_Porosity&, size_t)>;

const double* getDoublePr(const mxArray* arr)
{
    if (mxGetClassID(arr) != mxDOUBLE_CLASS)
        throw std::runtime_error("only double data handled");
    return mxGetPr(arr);
}
double* getDoublePr(mxArray* arr)
{
    if (mxGetClassID(arr) != mxDOUBLE_CLASS)
        throw std::runtime_error("only double data handled");
    return mxGetPr(arr);
}

class CommonArguments {
  public:
    CommonArguments(const mxArray** prhs, int nrhs)
        : mineral_property(getDoublePr(prhs[0])),
          fluid_property(getDoublePr(prhs[1])),
          phi_vector(getDoublePr(prhs[2])),
          len(mxGetM(prhs[0]))
    {
        switch (nrhs) {
            case (8):
                // mxGetScalar will convert if necessary
                frequency = mxGetScalar(prhs[4]);
                angle = mxGetScalar(prhs[5]);
                per_inc_con = mxGetScalar(prhs[6]);
                per_inc_ani = mxGetScalar(prhs[7]);
                v = nullptr;
                alpha = nullptr;
                break;
            case (9):
                // getDoublePr only allows double data
                alpha = getDoublePr(prhs[3]);
                v = getDoublePr(prhs[4]);
                frequency = mxGetScalar(prhs[5]);
                angle = mxGetScalar(prhs[6]);
                per_inc_con = mxGetScalar(prhs[7]);
                per_inc_ani = mxGetScalar(prhs[8]);
                break;
            default:
                throw std::invalid_argument("Wrong number of arguments");
        }
    }

    const double *mineral_property, *fluid_property, *phi_vector;
    const double *alpha, *v;
    double frequency, angle, per_inc_con, per_inc_ani;
    size_t len;

    std::array<double, 4> fluid(size_t n) const
    {
        return {{fluid_property[n], fluid_property[n + len],
                 fluid_property[n + 2 * len], fluid_property[n + 3 * len]}};
    }

    std::array<double, 3> mineral(size_t n) const
    {
        return {{mineral_property[n], mineral_property[n + len],
                 mineral_property[n + 2 * len]}};
    }
};

void set_err(bool const condition, const std::string& msg)
{
    if (condition) throw std::invalid_argument(msg);
}

template <typename F, typename P, typename Index>
bool equivalent(F&& /*fun*/, P&&, Index)
{
    return true;
}

template <typename F, typename P, typename Index, typename... Indices>
bool equivalent(F&& fun, P&& ptr, Index idx1, Index idx2, Indices&&... rest)
{
    return (fun(ptr[idx1]) == fun(ptr[idx2])) &&
           equivalent(fun, std::forward<P>(ptr), idx2,
                      std::forward<Indices>(rest)...);
}

void check_arguments(int nlhs, int nrhs, const mxArray* prhs[])
{
    set_err(nrhs != 8 && nrhs != 9, "Invalid number of input arguments.");
    set_err(
        nlhs != 1,
        "Invalid number of output arguments, T-matrix only has single output.");
    set_err(mxGetN(prhs[0]) != 3,
            "Mineral property array should be size (N, 3)");
    set_err(mxGetN(prhs[1]) != 4, "Fluid property array should be size (N, 4)");

    size_t offset = 0;
    if (nrhs == 9) offset += 1;

    set_err(offset == 0 && SCALAR != mxGetNumberOfElements(prhs[3]),
            "Invalid scenario argument.");

    set_err(offset != 0 && (mxGetN(prhs[3]) == 1 || mxGetN(prhs[4]) == 1),
            "Alpha and V need to be passed as (N x [2, 3, 4]) sized "
            "arrays");

    set_err(SCALAR != mxGetNumberOfElements(prhs[4 + offset]),
            "Invalid frequency argument.");
    set_err(SCALAR != mxGetNumberOfElements(prhs[5 + offset]),
            "Invalid angle argument.");

    if (SCALAR !=
        mxGetNumberOfElements(prhs[6 + offset]))  // check per_inc_con length
        set_err(!equivalent(mxGetM, prhs, 0, 1, 3, 4, 7, 8),
                "Arrays need to have equivalent row sizes!");
}

std::tuple<TMatrix_Porosity, tmatrix_eval> create_scenario_type(
    const mxArray** prhs, int nrhs)
{
    CommonArguments args(prhs, nrhs);

    scenario scen;

    switch (int(mxGetScalar(prhs[3]))) {
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
            throw std::invalid_argument(
                "scen parameter must be in:\n  1 - "
                "DUAL_POR_MOSTLY_ROUNDED,\n  2 - "
                "DUAL_POR_LITTLE_ROUNDED,\n  3 - MIXED_PORES,\n  4 - "
                "FLAT_PORES_AND_CRACKS.\n");
    }

    TMatrix_Porosity t_matrix(scen, args.per_inc_con, args.per_inc_ani);
    auto eval = [args](const TMatrix_Porosity& t_matrix, size_t n) -> Result {
        return t_matrix.evaluate(args.mineral(n), args.fluid(n),
                                 args.phi_vector[n], args.frequency,
                                 args.angle);
    };

    return std::make_tuple(std::move(t_matrix), std::move(eval));
}

std::tuple<TMatrix_Porosity, tmatrix_eval> create_const_alpha(
    const mxArray** prhs, int nrhs)
{
    CommonArguments args(prhs, nrhs);
    TMatrix_Porosity t_matrix(args.per_inc_con, args.per_inc_ani,
                              mxGetN(prhs[3]));

    std::vector<double> alpha;
    std::vector<double> v;

    alpha.insert(alpha.end(), getDoublePr(prhs[3]),
                 getDoublePr(prhs[3]) + mxGetNumberOfElements(prhs[3]));
    v.insert(v.end(), getDoublePr(prhs[4]),
             getDoublePr(prhs[4]) + mxGetNumberOfElements(prhs[4]));

    t_matrix.setAlpha(alpha);
    t_matrix.setV(v);

    auto eval = [args, alpha](const TMatrix_Porosity& t_matrix,
                              size_t n) -> Result {
        return t_matrix.evaluate(args.mineral(n), args.fluid(n),
                                 args.phi_vector[n], args.frequency,
                                 args.angle);
    };

    return std::make_tuple(std::move(t_matrix), std::move(eval));
}

std::tuple<TMatrix_Porosity, tmatrix_eval> create_all_dynamic(
    const mxArray** prhs, int nrhs)
{
    CommonArguments args(prhs, nrhs);
    size_t alpha_len = mxGetN(prhs[3]);

    TMatrix_Porosity t_matrix(args.per_inc_con, args.per_inc_ani, alpha_len);
    const double* per_inc_con = getDoublePr(prhs[7]);
    const double* per_inc_ani = getDoublePr(prhs[8]);
    auto eval = [args, alpha_len, per_inc_con, per_inc_ani](
                    TMatrix_Porosity& t_matrix, size_t n) -> Result {
        assert(args.alpha != nullptr);
        std::vector<double> alpha;
        std::vector<double> v;

        for (size_t i = 0; i < alpha_len; i++) {
            alpha.push_back(args.alpha[n + args.len * i]);
            v.push_back(args.v[n + args.len * i]);
        }

        return t_matrix.evaluate(
            per_inc_con[n], per_inc_ani[n], alpha, v, args.mineral(n),
            args.fluid(n), args.phi_vector[n], args.frequency, args.angle);
    };

    return std::make_tuple(std::move(t_matrix), std::move(eval));
}

std::tuple<TMatrix_Porosity, tmatrix_eval> create_dynamic_alpha(
    const mxArray** prhs, int nrhs)
{
    CommonArguments args(prhs, nrhs);
    size_t alpha_len = mxGetN(prhs[3]);

    TMatrix_Porosity t_matrix(args.per_inc_con, args.per_inc_ani, alpha_len);

    auto eval = [args, alpha_len](TMatrix_Porosity& t_matrix,
                                  size_t n) -> Result {
        assert(args.alpha != nullptr);
        std::vector<double> alpha;
        std::vector<double> v;

        for (size_t i = 0; i < alpha_len; i++) {
            alpha.push_back(args.alpha[n + args.len * i]);
            v.push_back(args.v[n + args.len * i]);
        }

        return t_matrix.evaluate(alpha, v, args.mineral(n), args.fluid(n),
                                 args.phi_vector[n], args.frequency,
                                 args.angle);
    };

    return std::make_tuple(std::move(t_matrix), std::move(eval));
}

enum class evaluation_t
{
    SCENARIOS,
    STATIC,
    DYNAMIC,
    ALL_DYNAMIC
};

evaluation_t get_evaluation_type(const mxArray** prhs, int nrhs)
{
  if (nrhs == 8) return evaluation_t::SCENARIOS;

    if (mxGetM(prhs[3]) == 1) return evaluation_t::STATIC;
    if (mxGetNumberOfElements(prhs[8]) > 1)
        return evaluation_t::ALL_DYNAMIC;
    else
        return evaluation_t::DYNAMIC;
}

/*This function emulates the behaviour of the T_matrix_porosity
 * MATLAB function
 */
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    try {
        check_arguments(nlhs, nrhs, prhs);

        TMatrix_Porosity Tmatrix;
        tmatrix_eval fun;

        size_t len = mxGetM(prhs[0]);

        switch (get_evaluation_type(prhs, nrhs)) {
            case evaluation_t::SCENARIOS:
                std::tie(Tmatrix, fun) = create_scenario_type(prhs, nrhs);
                break;
            case evaluation_t::STATIC:
                std::tie(Tmatrix, fun) = create_const_alpha(prhs, nrhs);
                break;
            case evaluation_t::DYNAMIC:
                std::tie(Tmatrix, fun) = create_dynamic_alpha(prhs, nrhs);
                break;
            case evaluation_t::ALL_DYNAMIC:
                std::tie(Tmatrix, fun) = create_all_dynamic(prhs, nrhs);
                break;
        }

        /* Write data back to MATLAB object */

        plhs[0] = mxCreateStructMatrix(1, 1, N_OUT_FLD, fieldnames);
        double* sources[N_OUT_FLD];

        for (size_t i = 0; i < N_OUT_FLD; i++) {
            mxArray* m = mxCreateDoubleMatrix(len, 1, mxREAL);
            sources[i] = getDoublePr(m);
            mxSetField(plhs[0], 0, fieldnames[i], m);
        }
        /*
#if defined(_OPENMP)
#pragma omp parallel default(shared)
#pragma omp for schedule(runtime)
#endif
*/
        for (size_t n = 0; n < len; n++) {
            Result res = fun(Tmatrix, n);

            sources[0][n] = res.pVelocity;
            sources[1][n] = res.horzShearVelocity;
            sources[2][n] = res.vertShearVelocity;
            sources[3][n] = res.frequency;
            sources[4][n] = res.rhoEffective;
        }
    } catch (std::exception& e) {
        mexErrMsgIdAndTxt("MATLAB:TMatrix_porosity_mwrapper", e.what());
    }
}
