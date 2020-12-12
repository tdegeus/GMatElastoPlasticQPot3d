/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

*/

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN2D_ELASTIC_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN2D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

inline Elastic::Elastic(double K, double G) :
    GMatElastic::Cartesian3d::Elastic(K, G)
{
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticQPot3d

#endif
