
# required to specify the c++ standard
cmake_minimum_required(VERSION 3.0)

# required for install
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# project settings
# ----------------

# name
project(xElastoPlasticQPot3d)

# file that contains the version information
set(parse_version include/xElastoPlasticQPot3d/ElastoPlasticQPot3d.h)

# header files
set(headers
  include/xElastoPlasticQPot3d/ElastoPlasticQPot3d.h
  include/xElastoPlasticQPot3d/Cartesian3d.hpp
  include/xElastoPlasticQPot3d/Matrix.hpp
  include/xElastoPlasticQPot3d/Elastic.hpp
  include/xElastoPlasticQPot3d/Cusp.hpp
  include/xElastoPlasticQPot3d/Smooth.hpp
)

# automatically parse the version number
file(READ "${parse_version}" version)
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT3D_WORLD_VERSION[ \t]+([0-9]+)" _ "${version}")
set(world_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT3D_MAJOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(major_version "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+ELASTOPLASTICQPOT3D_MINOR_VERSION[ \t]+([0-9]+)" _ "${version}")
set(minor_version "${CMAKE_MATCH_1}")
set(ELASTOPLASTICQPOT3D_VERSION_NUMBER ${world_version}.${major_version}.${minor_version})

# paths
# -----

set(ELASTOPLASTICQPOT3D_ROOT_DIR    "${CMAKE_INSTALL_PREFIX}")
set(ELASTOPLASTICQPOT3D_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(INCLUDE_INSTALL_DIR           "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(CMAKEPACKAGE_INSTALL_DIR      "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(PKGCONFIG_INSTALL_DIR         "${CMAKE_INSTALL_DATADIR}/pkgconfig")
set(fcmake                        "ElastoPlasticQPot3dConfig.cmake")
set(fpkg                          "ElastoPlasticQPot3d.pc")

# options
# -------

# configure pkg-config (default: on)
option(PKGCONFIG "Build pkg-config ${fpkg} file" ON)

# disable pkg-config for native Windows builds
if(WIN32 OR CMAKE_HOST_SYSTEM_NAME MATCHES Windows)
  option(PKGCONFIG "Build pkg-config ${fpkg} file" OFF)
endif()

# C++ standard
# ------------

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 14)
  set(CMAKE_CXX_STANDARD 14)
endif()

# configure CMake
# ---------------

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ElastoPlasticQPot3dConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/ElastoPlasticQPot3dConfig.cmake
  PATH_VARS ELASTOPLASTICQPOT3D_INCLUDE_DIR ELASTOPLASTICQPOT3D_ROOT_DIR
  INSTALL_DESTINATION ${CMAKEPACKAGE_INSTALL_DIR}
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# install
# -------

# pkg-config
if(PKGCONFIG)
  configure_file(${fpkg}.in ${fpkg} @ONLY)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fpkg} DESTINATION ${PKGCONFIG_INSTALL_DIR})
endif()

# CMake
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${fcmake} DESTINATION ${CMAKEPACKAGE_INSTALL_DIR})

# header files
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${headers} DESTINATION ${INCLUDE_INSTALL_DIR})

# print information to screen
# ---------------------------

message(STATUS "")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "|")
message(STATUS "| Use 'make install' to install in")
message(STATUS "|   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "|")
message(STATUS "| To specify a custom directory call")
message(STATUS "|   cmake /path/to/${PROJECT_NAME} -DCMAKE_INSTALL_PREFIX=yourprefix")
message(STATUS "|")
message(STATUS "| For custom paths, add the following line to your '~/.bashrc'")
message(STATUS "|   export PKG_CONFIG_PATH=${CMAKE_INSTALL_PREFIX}/share/pkgconfig:$PKG_CONFIG_PATH")
message(STATUS "|")
message(STATUS "+---------------------------------------------------------------------------------")
message(STATUS "")
