/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

*/

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

inline double Elastic::K() const
{
    return m_K;
}

inline double Elastic::G() const
{
    return m_G;
}

template <class T>
inline void Elastic::stress(const Tensor2& Eps, T&& Sig) const
{
    auto I = Cartesian3d::I2();
    auto epsm = trace(Eps) / 3.0;
    auto Epsd = Eps - epsm * I;
    xt::noalias(Sig) = 3.0 * m_K * epsm * I + 2.0 * m_G * Epsd;
}

inline Tensor2 Elastic::Stress(const Tensor2& Eps) const
{
    Tensor2 Sig;
    this->stress(Eps, Sig);
    return Sig;
}

inline double Elastic::energy(const Tensor2& Eps) const
{
    auto I = Cartesian3d::I2();
    auto epsm = trace(Eps) / 3.0;
    auto Epsd = Eps - epsm * I;
    auto epsd = std::sqrt(0.5 * A2_ddot_B2(Epsd, Epsd));
    auto U = 3.0 * m_K * std::pow(epsm, 2.0);
    auto V = 2.0 * m_G * std::pow(epsd, 2.0);
    return U + V;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticQPot3d

#endif
