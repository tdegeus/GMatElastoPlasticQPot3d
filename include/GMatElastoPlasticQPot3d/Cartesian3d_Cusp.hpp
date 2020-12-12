/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot3d

*/

#ifndef GMATELASTOPLASTICQPOT3D_CARTESIAN2D_CUSP_HPP
#define GMATELASTOPLASTICQPOT3D_CARTESIAN2D_CUSP_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticQPot3d {
namespace Cartesian3d {

inline Cusp::Cusp(double K, double G, const xt::xtensor<double, 1>& epsy, bool init_elastic)
    : m_K(K), m_G(G)
{
    xt::xtensor<double, 1> y = xt::sort(epsy);

    if (init_elastic) {
        if (y(0) != -y(1)) {
            y = xt::concatenate(xt::xtuple(xt::xtensor<double, 1>({-y(0)}), y));
        }
    }

    GMATELASTOPLASTICQPOT3D_ASSERT(y.size() > 1);

    m_yield = QPot::Static(0.0, y);
}

inline double Cusp::K() const
{
    return m_K;
}

inline double Cusp::G() const
{
    return m_G;
}

inline xt::xtensor<double, 1> Cusp::epsy() const
{
    return m_yield.yield();
}

inline auto Cusp::getQPot() const
{
    return m_yield;
}

inline auto* Cusp::refQPot()
{
    return &m_yield;
}

inline size_t Cusp::currentIndex() const
{
    return m_yield.currentIndex();
}

inline double Cusp::currentYieldLeft() const
{
    return m_yield.currentYieldLeft();
}

inline double Cusp::currentYieldRight() const
{
    return m_yield.currentYieldRight();
}

inline double Cusp::epsp() const
{
    return 0.5 * (m_yield.currentYieldLeft() + m_yield.currentYieldRight());
}

inline double Cusp::energy() const
{
    namespace GT = GMatTensor::Cartesian3d::pointer;
    std::array<double, 9> Epsd;
    double epsm = GT::Hydrostatic_deviatoric(&m_Eps[0], &Epsd[0]);
    double epsd = std::sqrt(0.5 * GT::A2s_ddot_B2s(&Epsd[0], &Epsd[0]));

    double U = 3.0 * m_K * std::pow(epsm, 2.0);

    double eps_min = 0.5 * (m_yield.currentYieldRight() + m_yield.currentYieldLeft());
    double deps_y = 0.5 * (m_yield.currentYieldRight() - m_yield.currentYieldLeft());

    double V = 2.0 * m_G * (std::pow(epsd - eps_min, 2.0) - std::pow(deps_y, 2.0));

    return U + V;
}

inline bool Cusp::checkYieldBoundLeft(size_t n) const
{
    return m_yield.checkYieldBoundLeft(n);
}

inline bool Cusp::checkYieldBoundRight(size_t n) const
{
    return m_yield.checkYieldBoundRight(n);
}

template <class T>
inline void Cusp::setStrainPtr(const T* arg)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;
    std::copy(arg, arg + 9, m_Eps.begin());

    std::array<double, 9> Epsd;
    double epsm = GT::hydrostatic_deviatoric(&m_Eps[0], &Epsd[0]);
    double epsd = std::sqrt(0.5 * GT::A2s_ddot_B2s(&Epsd[0], &Epsd[0]));
    m_yield.setPosition(epsd);

    m_Sig[0] = m_Sig[4] = m_Sig[8] = 3.0 * m_K * epsm;

    if (epsd <= 0.0) {
        m_Sig[1] = m_Sig[2] = m_Sig[3] = m_Sig[5] = m_Sig[6] = m_Sig[7] = 0.0;
        return;
    }

    double eps_min = 0.5 * (m_yield.currentYieldRight() + m_yield.currentYieldLeft());

    double g = 2.0 * m_G * (1.0 - eps_min / epsd);
    m_Sig[0] += g * Epsd[0];
    m_Sig[1] = g * Epsd[1];
    m_Sig[2] = g * Epsd[2];
    m_Sig[3] = g * Epsd[3];
    m_Sig[4] += g * Epsd[4];
    m_Sig[5] = g * Epsd[5];
    m_Sig[6] = g * Epsd[6];
    m_Sig[7] = g * Epsd[7];
    m_Sig[8] += g * Epsd[8];
}

template <class T>
inline void Cusp::strainPtr(T* ret) const
{
    std::copy(m_Eps.begin(), m_Eps.end(), ret);
}

template <class T>
inline void Cusp::stressPtr(T* ret) const
{
    std::copy(m_Sig.begin(), m_Sig.end(), ret);
}

template <class T>
inline void Cusp::tangentPtr(T* ret) const
{
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto C = m_K * II + 2.0 * m_G * I4d;
    std::copy(C.cbegin(), C.cend(), ret);
}

template <class T>
inline void Cusp::setStrain(const T& arg)
{
    GMATELASTOPLASTICQPOT3D_ASSERT(xt::has_shape(arg, {3, 3}));
    return this->setStrainPtr(arg.data());
}

template <class T>
inline void Cusp::strain(T& ret) const
{
    GMATELASTOPLASTICQPOT3D_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->strainPtr(ret.data());
}

template <class T>
inline void Cusp::stress(T& ret) const
{
    GMATELASTOPLASTICQPOT3D_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->stressPtr(ret.data());
}

template <class T>
inline void Cusp::tangent(T& ret) const
{
    GMATELASTOPLASTICQPOT3D_ASSERT(xt::has_shape(ret, {3, 3, 3, 3}));
    return this->tangentPtr(ret.data());
}

inline xt::xtensor<double, 2> Cusp::Strain() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->strainPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 2> Cusp::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->stressPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 4> Cusp::Tangent() const
{
    xt::xtensor<double, 4> ret = xt::empty<double>({3, 3, 3, 3});
    this->tangentPtr(ret.data());
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticQPot3d

#endif
