/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

================================================================================================= */

#ifndef GMATELASTOPLASTICQPOT3D_H
#define GMATELASTOPLASTICQPOT3D_H

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
#include <xtensor/xnoalias.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

// -------------------------------------------------------------------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------------------------------------------------------------------

#define GMATELASTOPLASTICQPOT3D_WORLD_VERSION 0
#define GMATELASTOPLASTICQPOT3D_MAJOR_VERSION 0
#define GMATELASTOPLASTICQPOT3D_MINOR_VERSION 3

#define GMATELASTOPLASTICQPOT3D_VERSION_AT_LEAST(x,y,z) \
  (GMATELASTOPLASTICQPOT3D_WORLD_VERSION>x || (GMATELASTOPLASTICQPOT3D_WORLD_VERSION>=x && \
  (GMATELASTOPLASTICQPOT3D_MAJOR_VERSION>y || (GMATELASTOPLASTICQPOT3D_MAJOR_VERSION>=y && \
                                               GMATELASTOPLASTICQPOT3D_MINOR_VERSION>=z))))

#define GMATELASTOPLASTICQPOT3D_VERSION(x,y,z) \
  (GMATELASTOPLASTICQPOT3D_WORLD_VERSION==x && \
   GMATELASTOPLASTICQPOT3D_MAJOR_VERSION==y && \
   GMATELASTOPLASTICQPOT3D_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.h"
#include "Cartesian3d.hpp"
#include "Cartesian3d_Elastic.hpp"
#include "Cartesian3d_Cusp.hpp"
#include "Cartesian3d_Smooth.hpp"
#include "Cartesian3d_Matrix.hpp"

// -------------------------------------------------------------------------------------------------

#endif
