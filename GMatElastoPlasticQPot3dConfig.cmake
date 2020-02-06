# GMatElastoPlasticQPot3d cmake module
#
# This module sets the target:
#
#     GMatElastoPlasticQPot3d
#
# In addition, it sets the following variables:
#
#     GMatElastoPlasticQPot3d_FOUND - true if GMatElastoPlasticQPot3d found
#     GMatElastoPlasticQPot3d_VERSION - GMatElastoPlasticQPot3d's version
#     GMatElastoPlasticQPot3d_INCLUDE_DIRS - the directory containing GMatElastoPlasticQPot3d headers
#
# The following support targets are defined to simplify things:
#
#     GMatElastoPlasticQPot3d::compiler_warnings - enable compiler warnings
#     GMatElastoPlasticQPot3d::assert - enable GMatElastoPlasticQPot3d assertions
#     GMatElastoPlasticQPot3d::debug - enable all assertions (slow)

include(CMakeFindDependencyMacro)

# Define target "GMatElastoPlasticQPot3d"

if(NOT TARGET GMatElastoPlasticQPot3d)
    include("${CMAKE_CURRENT_LIST_DIR}/GMatElastoPlasticQPot3dTargets.cmake")
    get_target_property(
        GMatElastoPlasticQPot3d_INCLUDE_DIRS
        GMatElastoPlasticQPot3d
        INTERFACE_INCLUDE_DIRECTORIES)
endif()

# Find dependencies

find_dependency(xtensor)

# Define support target "GMatElastoPlasticQPot3d::compiler_warnings"

if(NOT TARGET GMatElastoPlasticQPot3d::compiler_warnings)
    add_library(GMatElastoPlasticQPot3d::compiler_warnings INTERFACE IMPORTED)
    if(MSVC)
        set_property(
            TARGET GMatElastoPlasticQPot3d::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            /W4)
    else()
        set_property(
            TARGET GMatElastoPlasticQPot3d::compiler_warnings
            PROPERTY INTERFACE_COMPILE_OPTIONS
            -Wall -Wextra -pedantic -Wno-unknown-pragmas)
    endif()
endif()

# Define support target "GMatElastoPlasticQPot3d::assert"

if(NOT TARGET GMatElastoPlasticQPot3d::assert)
    add_library(GMatElastoPlasticQPot3d::assert INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticQPot3d::assert
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        GMATELASTIC_ENABLE_ASSERT)
endif()

# Define support target "GMatElastoPlasticQPot3d::debug"

if(NOT TARGET GMatElastoPlasticQPot3d::debug)
    add_library(GMatElastoPlasticQPot3d::debug INTERFACE IMPORTED)
    set_property(
        TARGET GMatElastoPlasticQPot3d::debug
        PROPERTY INTERFACE_COMPILE_DEFINITIONS
        XTENSOR_ENABLE_ASSERT GMATELASTIC_ENABLE_ASSERT)
endif()
