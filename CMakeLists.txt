cmake_minimum_required (VERSION 3.15 FATAL_ERROR)

project (Tmatrix)

# Require c++14
set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

option(USE_CONAN "Use Conan (https://conan.io) for package management" ON)

if(USE_CONAN)
    if(UNIX AND EXISTS "/etc/redhat-release")
        if(NOT DEFINED ENV{CONAN_CACERT_PATH})
            set(ENV{CONAN_CACERT_PATH} "/etc/pki/tls/cert.pem")
        elseif(NOT DEFINED ENV{REQUESTS_CA_BUNDLE})
            set(ENV{REQUESTS_CA_BUNDLE} "/etc/pki/tls/cert.pem")
        endif()
    endif()
endif()

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${CONAN_INCLUDE_DIRS_EIGEN})
add_definitions(-DEIGEN_NO_DEBUG -DEIGEN_MPL2_ONLY)
find_package(pybind11 REQUIRED)

add_subdirectory(src/tmatrix)

