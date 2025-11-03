# tmatrix

[![SCM Compliance](https://scm-compliance-api.radix.equinor.com/repos/equinor/tmatrix/badge)](https://developer.equinor.com/governance/scm-policy/)
[![PyPI version](https://badge.fury.io/py/tmatrix.svg)](https://badge.fury.io/py/tmatrix)

## Installation

For most users, installing from PyPI is the preferred way:

```bash
pip install tmatrix
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
