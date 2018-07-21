/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot3d

================================================================================================= */

#ifndef ELASTOPLASTICQPOT3D_ELASTIC_HPP
#define ELASTOPLASTICQPOT3D_ELASTIC_HPP

// -------------------------------------------------------------------------------------------------

#include "ElastoPlasticQPot3d.h"

// =================================================================================================

namespace ElastoPlasticQPot3d {

// ------------------------------------------ constructor ------------------------------------------

inline Elastic::Elastic(double K, double G) : m_kappa(K), m_mu(G)
{
}

// ------------------------------------------ parameters -------------------------------------------

inline double Elastic::kappa() const
{
  return m_kappa;
}

// ------------------------------------------ parameters -------------------------------------------

inline double Elastic::mu() const
{
  return m_mu;
}

// ---------------------------------- equivalent deviator strain -----------------------------------

inline double Elastic::epsd(const T2s &Eps) const
{
  T2s Epsd = Eps - Eps.trace()/3. * T2d::I();

  return std::sqrt(.5*Epsd.ddot(Epsd));
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Elastic::epsp(const T2s &Eps) const
{
  UNUSED(Eps);

  return 0.0;
}

// ----------------------------------- equivalent plastic strain -----------------------------------

inline double Elastic::epsp(double epsd) const
{
  UNUSED(epsd);

  return 0.0;
}

// ----------------------------------------- yield stress ------------------------------------------

inline double Elastic::epsy(size_t i) const
{
  UNUSED(i);

  return std::numeric_limits<double>::infinity();
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Elastic::find(const T2s &Eps) const
{
  UNUSED(Eps);

  return 0;
}

// ------------------------------------- find potential index --------------------------------------

inline size_t Elastic::find(double epsd) const
{
  UNUSED(epsd);

  return 0;
}

// -------------------------------------------- stress ---------------------------------------------

inline T2s Elastic::Sig(const T2s &Eps) const
{
  // decompose strain
  T2d    I     = T2d::I();
  double treps = Eps.trace();
  T2s    Epsd  = Eps - (treps/3.) * I;

  // return stress tensor
  return m_kappa * treps * I + 2. * m_mu * Epsd;
}

// -------------------------------------------- energy ---------------------------------------------

inline double Elastic::energy(const T2s &Eps) const
{
  // decompose strain
  double treps = Eps.trace();
  T2s    Epsd  = Eps - (treps/3.) * T2d::I();
  double epsd  = std::sqrt(.5*Epsd.ddot(Epsd));

  // energy components
  double U = 0.5 * m_kappa * std::pow(treps,2.);
  double V = 2.0 * m_mu    * std::pow(epsd ,2.);

  // return total energy
  return U + V;
}

// =================================================================================================

} // namespace ...

// =================================================================================================

#endif
