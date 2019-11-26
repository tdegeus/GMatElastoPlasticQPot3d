/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

================================================================================================= */

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN2D_CUSP_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN2D_CUSP_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Cusp::Cusp(double K, double G, const xt::xtensor<double,1>& epsy, bool init_elastic) :
  m_K(K), m_G(G)
{
  m_epsy = xt::sort(epsy);

  if (init_elastic)
    if (m_epsy(0) != -m_epsy(1))
      m_epsy = xt::concatenate(xt::xtuple(xt::xtensor<double,1>({-m_epsy(0)}), m_epsy));

  GMATELASTOPLASTICQPOT3D_ASSERT(m_epsy.size() > 1);
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::epsp(const Tensor2& Eps) const
{
  return this->epsp(Cartesian3d::Epsd(Eps));
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::epsp(double epsd) const
{
  size_t i = this->find(epsd);
  return 0.5 * ( m_epsy(i+1) + m_epsy(i) );
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::epsy(size_t i) const
{
  return m_epsy(i);
}

// -------------------------------------------------------------------------------------------------

inline size_t Cusp::find(const Tensor2& Eps) const
{
  return this->find(Cartesian3d::Epsd(Eps));
}

// -------------------------------------------------------------------------------------------------

inline size_t Cusp::find(double epsd) const
{
  GMATELASTOPLASTICQPOT3D_ASSERT(epsd > m_epsy(0) && epsd < m_epsy(m_epsy.size()-1));

  return std::lower_bound(m_epsy.begin(), m_epsy.end(), epsd) - m_epsy.begin() - 1;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void Cusp::stress(const Tensor2& Eps, T&& Sig) const
{
  // decompose strain: hydrostatic part, deviatoric part
  auto I = Cartesian3d::I2();
  auto epsm = trace(Eps) / 3.0;
  auto Epsd = Eps - epsm * I;
  auto epsd = std::sqrt(0.5 * A2_ddot_B2(Epsd, Epsd));

  // no deviatoric strain -> only hydrostatic stress
  if (epsd <= 0.) {
    xt::noalias(Sig) = 3.0 * m_K * epsm * I;
    return;
  }

  // read current yield strains
  size_t i = this->find(epsd);
  double eps_min = 0.5 * (m_epsy(i+1) + m_epsy(i));

  // return stress tensor
  xt::noalias(Sig) = 3.0 * m_K * epsm * I + 2.0 * m_G * (1.0 - eps_min/epsd) * Epsd;
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Cusp::Stress(const Tensor2& Eps) const
{
  Tensor2 Sig;
  this->stress(Eps, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

inline double Cusp::energy(const Tensor2& Eps) const
{
  // decompose strain: hydrostatic part, deviatoric part
  auto I = Cartesian3d::I2();
  auto epsm = trace(Eps) / 3.0;
  auto Epsd = Eps - epsm * I;
  auto epsd = std::sqrt(0.5 * A2_ddot_B2(Epsd, Epsd));

  // hydrostatic part of the energy
  double U = 3.0 * m_K * std::pow(epsm, 2.0);

  // read current yield strain
  size_t i = this->find(epsd);
  double eps_min = 0.5 * (m_epsy(i+1) + m_epsy(i));
  double deps_y = 0.5 * (m_epsy(i+1) - m_epsy(i));

  // deviatoric part of the energy
  double V = 2.0 * m_G * (std::pow(epsd - eps_min, 2.0) - std::pow(deps_y, 2.0));

  // return total energy
  return U + V;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
