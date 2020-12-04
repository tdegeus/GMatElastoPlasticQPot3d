
cmake_minimum_required(VERSION 3.0)

# Basic settings
# ==============

project(GMatElastoPlasticQPot3d)

option(BUILD_TESTS "${PROJECT_NAME} Build tests" OFF)

# Version
# =======

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/GMatElastoPlasticQPot3d/config.h" _GMAT_ver
     REGEX "#define GMATELASTOPLASTICQPOT3D_VERSION_(MAJOR|MINOR|PATCH)")

foreach(ver ${_GMAT_ver})
    if(ver MATCHES "#define GMATELASTOPLASTICQPOT3D_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(_GMAT_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()

set(_GMAT_VERSION ${_GMAT_VERSION_MAJOR}.${_GMAT_VERSION_MINOR}.${_GMAT_VERSION_PATCH})

message(STATUS "Building GMatElastoPlasticQPot3d v${_GMAT_VERSION}")

# Set target
# ==========

find_package(xtensor REQUIRED)
find_package(GMatTensor REQUIRED)
find_package(GMatElastic REQUIRED)
find_package(QPot REQUIRED)

add_library(GMatElastoPlasticQPot3d INTERFACE)

target_include_directories(GMatElastoPlasticQPot3d INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(GMatElastoPlasticQPot3d INTERFACE xtensor QPot GMatTensor GMatElastic)

# Installation
# ============

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/" DESTINATION include)

install(TARGETS GMatElastoPlasticQPot3d EXPORT GMatElastoPlasticQPot3d-targets)

install(
    EXPORT GMatElastoPlasticQPot3d-targets
    FILE GMatElastoPlasticQPot3dTargets.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlasticQPot3d")

set(_GMAT ${CMAKE_SIZEOF_VOID_P})
unset(CMAKE_SIZEOF_VOID_P)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticQPot3dConfigVersion.cmake"
    VERSION ${_GMAT_VERSION}
    COMPATIBILITY AnyNewerVersion)

set(CMAKE_SIZEOF_VOID_P ${_GMAT})

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticQPot3dConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticQPot3dConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/GMatElastoPlasticQPot3d")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GMatElastoPlasticQPot3d.pc.in"
               "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticQPot3d.pc" @ONLY)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/GMatElastoPlasticQPot3d.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig/")

# Add builds
# ==========

include("GMatElastoPlasticQPot3dConfig.cmake")

set(CTEST_TEST_TARGET_ALIAS RUN_ALL_TESTS)
include(CTest)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
