# tmatrix

[![SCM Compliance](https://scm-compliance-api.radix.equinor.com/repos/equinor/tmatrix/badge)](https://developer.equinor.com/governance/scm-policy/)
[![PyPI version](https://badge.fury.io/py/tmatrix.svg)](https://badge.fury.io/py/tmatrix)

## Installation

For most users, installing from PyPI is the preferred way:

```bash
pip install tmatrix
```

Or using `uv`:

```bash
uv add tmatrix
```

For developers, the project can be compiled with `cmake`:

```bash
cd tmatrix
mkdir build && cd build
cmake ..
make
```

All objects are placed in the build subdirectory.

Note that enabling parallel processing incurs some overhead, and should only be
enabled for large jobs (e.g. 10.000+ sequential calls).

Under Windows use, find your desired Windows CMake [generator](https://cmake.org/cmake/help/v3.4/manual/cmake-generators.7.html#visual-studio-generators), ie:

```bash
cd tmatrix
mkdir build
cd build
cmake .. -G "Visual Studio 14 2015 Win64"
cmake --build . --target ALL_BUILD --config Release
```

## Usage

The package exposes two functions

- `tmatrix_porosity`
- `tmatrix_porosity_noscenario`

### TMatrix porosity

```python
from tmatrix import tmatrix_porosity


# Dimension of the output array
dim = 21

# Output result is stored in `out_np`
out_np = np.zeros((dim, 4))

# Mineral properties. Contains mineral bulk modulus [Pa], shear modulus [Pa] and density [kg/m³]. Shape should be (N, 3).
mineral_property_np = np.tile(np.array([7.10e10, 3.20e10, 2.71e03]), (dim, 1))

# Mineral properties. Contains mineral bulk modulus [Pa], shear modulus [Pa] and density [kg/m³]. Shape should be (N, 3).  
fluid_property_np = np.tile(np.array([2.700e09, 1.005e03, 1.000e02, 1.000e02]), (dim, 1))

# Porosity values. Shape should be (N,).
phi_vector_np = np.linspace(0.15, 0.25, dim)

# Input scenario. Can be 1,2,3 or 4.
#   1: Dual porosity, mostly rounded pores
#   2: Dual porosity, little rounded pores
#   3: Mixed pores
#   4: Flat pores and cracks
in_scenario = 1

# Signal frequency [Hz]
frequency = 1000

# Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg]
angle_of_sym_plane = 90

# Fraction of inclusions that are connected
per_inc_con = 0.5

# Fraction of inclusions that are anisotropic
per_inc_any = 0.5

_ = tmatrix.tmatrix_porosity(
    out_np=out_np,
    dim=dim,
    mineral_property_np=mineral_property_np,
    fluid_property_np=fluid_property_np,
    phi_vector_np=phi_vector_np,
    in_scenario=in_scenario,
    frequency=frequency,
    angle_of_sym_plane=angle_of_sym_plane,
    per_inc_con=per_inc_con,
    per_inc_any=per_inc_any,
)

# Returns 0 if success, otherwise failure. Result will be stored in `out_np`, with shape (dim, 4). 
# Column values in order are:
#   Vp: Vertical P-wave velocity [m/s]
#   Vsv: Vertical polarity S-wave velocity [m/s]
#   Vsh: Horizontal polarity S-wave velocity [m/s]
#   Rhob [kg/m^3]
```

### TMatrix porosity noscenario

```python
from tmatrix import tmatrix_porosity_noscenario

# Dimension of the output array
out_N = 21

# Output result is stored in `out_np`
out_np = np.zeros((out_N, 4))

# Mineral properties. Contains mineral bulk modulus [Pa], shear modulus [Pa] and density [kg/m³]. Shape should be (N, 3).
mineral_property_np = np.tile(np.array([7.10e10, 3.20e10, 2.71e03]), (out_N, 1))

# Fluid properties. Contains fluid bulk modulus [Pa] and density [kg/m³], viscosity [cP] and permeability [mD]. Shape should be (N, 4).
fluid_property_np = np.tile(np.array([2.700e09, 1.005e03, 1.000e02, 1.000e02]), (out_N, 1))

# Porosity values. Shape should be (N,).
phi_vector_np = np.linspace(0.15, 0.25, out_N)

# Aspect ratio values. Shape should be (N,) where N is the number of aspect ratio values
alpha_np = np.tile(np.array([0.9, 0.1]), (out_N, 1))

# Number of aspect ratio values per sample
alpha_size_np = np.full((out_N,), 2, dtype=int)

# Length of alpha array
alpha_N = 21

# Fraction of porosity with given aspect ratio
v_np = np.tile(np.array([0.9, 0.1]), (out_N, 1))

# Signal frequency [Hz]
frequency = 1000

# Angle of symmetry plane (0 = HTI, 90 = VTI medium) [deg]
angle = 90

# Fraction of inclusions that are connected
inc_con_np = np.array([0,5])

# Fraction of inclusions that are anisotropic
inc_ani_np = np.array([0,5])

# Length of `inc_con_np` and `inc_ani_np` 
inc_con_N = 1

tmatrix.tmatrix_porosity_noscenario(
    out_np=out_np,
    out_N=out_N,
    mineral_property_np=mineral_property_np,
    fluid_property_np=fluid_property_np,
    phi_vector_np=phi_vector_np,
    alpha_np=alpha_np,
    v_np=v_np,
    alpha_size_np=alpha_size_np,
    alpha_N=alpha_N,
    frequency=frequency,
    angle=angle,
    inc_con_np=inc_con_np,
    inc_ani_np=inc_ani_np,
    inc_con_N=inc_con_N,
)
# Returns None. Result will be stored in `out_np`. Output array has shape (out_N, 4).
# Column values in order are:
#   Vp: Vertical P-wave velocity [m/s]
#   Vsv: Vertical polarity S-wave velocity [m/s]
#   Vsh: Horizontal polarity S-wave velocity [m/s]
#   Rhob [kg/m^3]

```

## Literature

The theory can be found in the papers and in the references therein:

1. Agersborg, R., Jakobsen, M., Ruud, B.O. and Johansen, T. A. 2007.
Effects of pore fluid pressure on the seismic response of a fractured carbonate reservoir.
Stud. Geophys. Geod., 51, 89-118.
[Link](dx.doi.org/10.1007/s11200-007-0005-8)

2. Agersborg, R., Johansen, T. A. and Ruud, B.O. 2008.
Modelling reflection signatures of pore fluids and dual porosity in carbonate reservoirs.
Journal of Seismic Exploration, 17(1), 63-83.

3. Agersborg, R., Johansen, T. A., Jakobsen, M., Sothcott, J. and Best, A. 2008.
Effect of fluids and dual-pores systems on pressure-dependent velocities and attenuation in carbonates,
Geophysics, 73, No. 5, N35-N47.
[Link](dx.doi.org/10.1190/1.2969774)

4. Agersborg, R., Johansen, T. A., and Jakobsen, M. 2009.
Velocity variations in carbonate rocks due to dual porosity and wave-induced fluid flow.
Geophysical Prospecting, 57, 81-98.
[Link](dx.doi.org/10.1111/j.1365-2478.2008.00733.x)

All of the papers and a extended explanations of the involved equations
can be found in  Agersborg (2007), phd thesis:
[Link](https://bora.uib.no/handle/1956/2422)
