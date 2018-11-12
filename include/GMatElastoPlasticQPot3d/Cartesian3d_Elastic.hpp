/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

================================================================================================= */

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN3D_ELASTIC_HPP

// -------------------------------------------------------------------------------------------------

#include "Cartesian3d.h"

// =================================================================================================

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Elastic::Elastic(double K, double G) : m_kappa(K), m_mu(G)
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsd(const T2s &Eps) const
{
  auto Epsd = Eps - trace(Eps)/3. * eye();

  return std::sqrt(.5*ddot(Epsd,Epsd));
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsp(const T2s &) const
{
  return 0.0;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsp(double) const
{
  return 0.0;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsy(size_t) const
{
  return std::numeric_limits<double>::infinity();
}

// -------------------------------------------------------------------------------------------------

inline size_t Elastic::find(const T2s &) const
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

inline size_t Elastic::find(double) const
{
  return 0;
}

// -------------------------------------------------------------------------------------------------

inline T2s Elastic::Sig(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I     = eye();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;

  // return stress tensor
  return m_kappa * treps * I + 2.0 * m_mu * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I     = eye();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;
  auto epsd  = std::sqrt(.5*ddot(Epsd,Epsd));

  // hydrostatic part of the energy
  double U = 0.5 * m_kappa * std::pow(treps,2.);
  // deviatoric part of the energy
  double V = 2.0 * m_mu    * std::pow(epsd ,2.);

  // return total energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif
