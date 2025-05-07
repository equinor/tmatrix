set(SRC "${CMAKE_CURRENT_LIST_DIR}/src")

set (sources
    ${SRC}/calc_C_eff.cpp
    ${SRC}/calc_C_eff_visco.cpp
    ${SRC}/calc_isolatedpart.cpp
    ${SRC}/calc_Kd.cpp
    ${SRC}/calc_Kd_uuvv.cpp
    ${SRC}/calc_t.cpp
    ${SRC}/calc_td.cpp
    ${SRC}/calc_Theta.cpp
    ${SRC}/calc_X.cpp
    ${SRC}/calc_Z.cpp
    ${SRC}/Gtensor.cpp
    ${SRC}/iso_av_all.cpp
    ${SRC}/iso_average.cpp
    ${SRC}/TMatrix_porosity.cpp
    ${SRC}/velocity_vti_angles.cpp
)

pybind11_add_module(_tmatrix "${SRC}/python.cpp" ${sources})

target_link_libraries(_tmatrix PRIVATE Eigen3::Eigen)

# Installation directory is set in pyproject.toml
install(TARGETS _tmatrix LIBRARY DESTINATION .)
