/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

*/

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN3D_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

template <class T, class U>
inline void epsd(const T& A, U& B)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, B);
    B *= std::sqrt(0.5);
}

template <class T>
inline auto Epsd(const T& A)
{
    return xt::eval(std::sqrt(0.5) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

template <class T, class U>
inline void sigd(const T& A, U& B)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, B);
    B *= std::sqrt(2.0);
}

template <class T>
inline auto Sigd(const T& A)
{
    return xt::eval(std::sqrt(2.0) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticQPot3d

#endif
