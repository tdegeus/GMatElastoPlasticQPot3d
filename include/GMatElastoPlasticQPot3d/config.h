/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

================================================================================================= */

#ifndef GMATELASTOPLASTICQPOT3D_CONFIG_H
#define GMATELASTOPLASTICQPOT3D_CONFIG_H

// -------------------------------------------------------------------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

#ifndef NDEBUG
#define GMATELASTOPLASTICQPOT3D_ENABLE_ASSERT
#endif

#ifdef GMATELASTOPLASTICQPOT3D_ENABLE_ASSERT
#define GMATELASTOPLASTICQPOT3D_ASSERT(expr) GMATELASTOPLASTICQPOT3D_ASSERT_IMPL(expr, __FILE__, __LINE__)
#define GMATELASTOPLASTICQPOT3D_ASSERT_IMPL(expr, file, line)                                                               \
    if (!(expr))                                                                                                          \
    {                                                                                                                     \
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": assertion failed (" #expr ") \n\t"); \
    }
#else
#define GMATELASTOPLASTICQPOT3D_ASSERT(expr)
#endif

// -------------------------------------------------------------------------------------------------

#define GMATELASTOPLASTICQPOT3D_VERSION_MAJOR 0
#define GMATELASTOPLASTICQPOT3D_VERSION_MINOR 0
#define GMATELASTOPLASTICQPOT3D_VERSION_PATCH 4

#define GMATELASTOPLASTICQPOT3D_VERSION_AT_LEAST(x,y,z) \
  (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR > x || (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR >= x && \
  (GMATELASTOPLASTICQPOT3D_VERSION_MINOR > y || (GMATELASTOPLASTICQPOT3D_VERSION_MINOR >= y && \
                                                 GMATELASTOPLASTICQPOT3D_VERSION_PATCH >= z))))

#define GMATELASTOPLASTICQPOT3D_VERSION(x,y,z) \
  (GMATELASTOPLASTICQPOT3D_VERSION_MAJOR == x && \
   GMATELASTOPLASTICQPOT3D_VERSION_MINOR == y && \
   GMATELASTOPLASTICQPOT3D_VERSION_PATCH == z)

// -------------------------------------------------------------------------------------------------

#endif
