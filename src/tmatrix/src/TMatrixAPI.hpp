#pragma once

#include <array>
#include <vector>

enum class scenario
{
    DUAL_POR_MOSTLY_ROUNDED = 0,
    DUAL_POR_LITTLE_ROUNDED = 1,
    MIXED_PORES = 2,
    FLAT_PORES_AND_CRACKS = 3,
    LAST_SCENARIO = FLAT_PORES_AND_CRACKS
    //    NO_SCENARIO = -1
};

enum class control1
{
    ISOTROPIC = 0,
    MIXED_CASE = 1,
    ANISOTROPIC = 2
};

enum class control2
{
    ISOLATED = 0,
    BOTH = 1,
    CONNECTED = 2
};

struct Result
{
    double pVelocity;
    double horzShearVelocity;
    double vertShearVelocity;
    double frequency;
    double rhoEffective;
    double gamma;
    double delta;
    double epsilon;
};

class TMatrix_Porosity {
    std::vector<double> alpha;
    std::vector<double> v;
    std::vector<double> tau;
    control1 ctrl1;
    control2 ctrl2;
    double per_inc_con;
    double per_inc_ani;

    void updateControl(void);

  public:
    // create empty object
    TMatrix_Porosity(){};
    TMatrix_Porosity(scenario scen, double per_inc_con, double per_inc_ani);
    TMatrix_Porosity(double per_inc_con, double per_inc_ani, size_t len);

    void setAlpha(std::vector<double> new_alpha);

    void setV(const std::vector<double> new_v);

    Result evaluate(const std::array<double, 3>& mineralProperty,
                    const std::array<double, 4>& fluidProperty,
                    double const porosity, double const frequency,
                    double const angle) const;

    template <typename... Ts>
    Result evaluate(double con, double ani, Ts&&... args)
    {
        per_inc_con = con;
        per_inc_ani = ani;
        updateControl();

        return evaluate(std::forward<Ts>(args)...);
    }

    template <typename... Ts>
    Result evaluate(const std::vector<double>& alpha_,
                    const std::vector<double>& v_, Ts&&... args)
    {
        setAlpha(alpha_);
        setV(v_);
        return evaluate(std::forward<Ts>(args)...);
    }
};
