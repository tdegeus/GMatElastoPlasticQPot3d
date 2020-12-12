# GMatElastoPlasticQPot3d cmake module
#
# This module sets the target:
#
#     GMatElastoPlasticQPot3d
#
# In addition, it sets the following variables:
#
#     GMatElastoPlasticQPot3d_FOUND - true if the library is found
#     GMatElastoPlasticQPot3d_VERSION - the library's version
#     GMatElastoPlasticQPot3d_INCLUDE_DIRS - directory containing the library's headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastoPlasticQPot3d::compiler_warnings - enable compiler warnings
#     GMatElastoPlasticQPot3d::assert - enable library assertions
#     GMatElastoPlasticQPot3d::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastoPlasticQPot3d"

if(NOT TARGET GMatElastoPlasticQPot3d)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElastoPlasticQPot3dTargets.cmake")
endif()

# Define "GMatElastoPlasticQPot3d_INCLUDE_DIRS"

get_target_property(
    GMatElastoPlasticQPot3d_INCLUDE_DIRS
    GMatElastoPlasticQPot3d
    INTERFACE_INCLUDE_DIRECTORIES)

# Find dependencies

find_dependency(xtensor)
find_dependency(QPot)
find_dependency(GMatTensor)
find_dependency(GMatElastic)

# Define support target "GMatElastoPlasticQPot3d::compiler_warnings"

if(NOT TARGET GMatElastoPlasticQPot3d::compiler_warnings)
    add_library(GMatElastoPlasticQPot3d::compiler_warnings INTERFACE IMPORTED)
    target_link_libraries(GMatElastoPlasticQPot3d::compiler_warnings INTERFACE
        GMatTensor::compiler_warnings)
endif()

# Define support target "GMatElastoPlasticQPot3d::assert"

if(NOT TARGET GMatElastoPlasticQPot3d::assert)
    add_library(GMatElastoPlasticQPot3d::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticQPot3d::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTOPLASTICQPOT3D_ENABLE_ASSERT
        GMATTENSOR_ENABLE_ASSERT
        QPOT_ENABLE_ASSERT)
endif()

# Define support target "GMatElastoPlasticQPot3d::debug"

if(NOT TARGET GMatElastoPlasticQPot3d::debug)
    add_library(GMatElastoPlasticQPot3d::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticQPot3d::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT
        GMATELASTOPLASTICQPOT3D_ENABLE_ASSERT
        GMATTENSOR_ENABLE_ASSERT
        QPOT_ENABLE_ASSERT)
endif()
