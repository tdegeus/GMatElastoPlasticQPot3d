/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

================================================================================================= */

#ifndef XGMATELASTOPLASTICQPOT3D_CARTESIAN3D_SMOOTH_HPP
#define XGMATELASTOPLASTICQPOT3D_CARTESIAN3D_SMOOTH_HPP

// -------------------------------------------------------------------------------------------------

#include "GMatElastoPlasticQPot3d.h"

// =================================================================================================

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Smooth::Smooth(double kappa, double mu, const xt::xtensor<double,1> &epsy, bool init_elastic) :
  m_kappa(kappa), m_mu(mu)
{
  // copy sorted yield strains
  m_epsy = xt::sort(epsy);

  // extra yield strain, to force an initial elastic response
  if ( init_elastic )
    if ( m_epsy(0) != -m_epsy(1) )
      m_epsy = xt::concatenate(xt::xtuple(xt::xtensor<double,1>({-m_epsy(0)}), m_epsy));

  // check the number of yield strains
  if ( m_epsy.size() < 2 )
    throw std::runtime_error("Specify at least two yield strains 'epsy'");
}

// -------------------------------------------------------------------------------------------------

inline double Smooth::epsd(const T2s &Eps) const
{
  auto Epsd = Eps - trace(Eps)/3. * eye();

  return std::sqrt(.5*ddot(Epsd,Epsd));
}

// -------------------------------------------------------------------------------------------------

inline double Smooth::epsp(const T2s &Eps) const
{
  return epsp(epsd(Eps));
}

// -------------------------------------------------------------------------------------------------

inline double Smooth::epsp(double epsd) const
{
  size_t i = find(epsd);

  return ( m_epsy(i+1) + m_epsy(i) ) / 2.;
}

// -------------------------------------------------------------------------------------------------

inline double Smooth::epsy(size_t i) const
{
  return m_epsy(i);
}

// -------------------------------------------------------------------------------------------------

inline size_t Smooth::find(const T2s &Eps) const
{
  return find(epsd(Eps));
}

// -------------------------------------------------------------------------------------------------

inline size_t Smooth::find(double epsd) const
{
  if ( epsd <= m_epsy(0) or epsd >= m_epsy(m_epsy.size()-1) )
    throw std::runtime_error("Insufficient 'epsy'");

  return std::lower_bound(m_epsy.begin(), m_epsy.end(), epsd) - m_epsy.begin() - 1;
}

// -------------------------------------------------------------------------------------------------

inline T2s Smooth::Sig(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I     = eye();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;
  auto epsd  = std::sqrt(.5*ddot(Epsd,Epsd));

  // no deviatoric strain -> only hydrostatic stress
  if ( epsd <= 0. ) return m_kappa * treps * I;

  // read current yield strains
  size_t i       = find(epsd);
  double eps_min = ( m_epsy(i+1) + m_epsy(i) ) / 2.;
  double deps_y  = ( m_epsy(i+1) - m_epsy(i) ) / 2.;

  // return stress tensor
  return m_kappa*treps*I + (2.*m_mu/epsd)*(deps_y/M_PI)*sin(M_PI/deps_y*(epsd-eps_min))*Epsd;
}

// -------------------------------------------------------------------------------------------------

inline double Smooth::energy(const T2s &Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  T2s  I     = eye();
  auto treps = trace(Eps);
  auto Epsd  = Eps - treps/3. * I;
  auto epsd  = std::sqrt(.5*ddot(Epsd,Epsd));

  // hydrostatic part of the energy
  double U = 0.5 * m_kappa * std::pow(treps,2.);

  // read current yield strain
  size_t i       = find(epsd);
  double eps_min = ( m_epsy(i+1) + m_epsy(i) ) / 2.;
  double deps_y  = ( m_epsy(i+1) - m_epsy(i) ) / 2.;

  // deviatoric part of the energy
  double V = -4.0 * m_mu * std::pow(deps_y/M_PI,2.) * ( 1. + cos( M_PI/deps_y * (epsd-eps_min) ) );

  // return total energy
  return U + V;
}

// =================================================================================================

}} // namespace ...

// =================================================================================================

#endif
