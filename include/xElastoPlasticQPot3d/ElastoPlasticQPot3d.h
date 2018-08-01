/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/ElastoPlasticQPot3d

================================================================================================= */

#ifndef XELASTOPLASTICQPOT3D_H
#define XELASTOPLASTICQPOT3D_H

// --------------------------------------- include libraries ---------------------------------------

// use "M_PI" from "math.h"
#define _USE_MATH_DEFINES

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>

// ---------------------------------------- dummy operation ----------------------------------------

// dummy operation that can be use to suppress the "unused parameter" warnings
#define UNUSED(p) ( (void)(p) )

// -------------------------------------- version information --------------------------------------

#define ELASTOPLASTICQPOT3D_WORLD_VERSION 0
#define ELASTOPLASTICQPOT3D_MAJOR_VERSION 0
#define ELASTOPLASTICQPOT3D_MINOR_VERSION 2

#define ELASTOPLASTICQPOT3D_VERSION_AT_LEAST(x,y,z) \
  (ELASTOPLASTICQPOT3D_WORLD_VERSION>x || (ELASTOPLASTICQPOT3D_WORLD_VERSION>=x && \
  (ELASTOPLASTICQPOT3D_MAJOR_VERSION>y || (ELASTOPLASTICQPOT3D_MAJOR_VERSION>=y && \
                                           ELASTOPLASTICQPOT3D_MINOR_VERSION>=z))))

#define ELASTOPLASTICQPOT3D_VERSION(x,y,z) \
  (ELASTOPLASTICQPOT3D_WORLD_VERSION==x && \
   ELASTOPLASTICQPOT3D_MAJOR_VERSION==y && \
   ELASTOPLASTICQPOT3D_MINOR_VERSION==z)

// ====================================== ElastoPlasticQPot3d ======================================

namespace xElastoPlasticQPot3d {

// --------------------------------------------- alias ---------------------------------------------

using T2s = xt::xtensor_fixed<double, xt::xshape<3,3>>;

// ---------------------------------------- tensor algebra -----------------------------------------

template<class T> inline double trace(const T &A);
template<class T> inline double ddot (const T &A, const T &B);

// -------------------------- equivalent stress/strain (Cartesian3d.cpp) ---------------------------

// mean
inline double sigm(const T2s &Sig);
inline double epsm(const T2s &Eps);

// equivalent deviator
inline double sigd(const T2s &Sig);
inline double epsd(const T2s &Eps);

// deviator
inline T2s Sigd(const T2s &Sig);
inline T2s Epsd(const T2s &Eps);

// ----------------------------- equivalent stress/strain (Matrix.cpp) -----------------------------

// no allocation
inline void sigm(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigm);
inline void epsm(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsm);
inline void sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,2> &a_sigd);
inline void epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsd);
inline void Sigd(const xt::xtensor<double,4> &a_Sig, xt::xtensor<double,4> &a_Sigd);
inline void Epsd(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Epsd);

// return allocated result
inline xt::xtensor<double,2> sigm(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,2> epsm(const xt::xtensor<double,4> &a_Eps);
inline xt::xtensor<double,2> sigd(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,2> epsd(const xt::xtensor<double,4> &a_Eps);
inline xt::xtensor<double,4> Sigd(const xt::xtensor<double,4> &a_Sig);
inline xt::xtensor<double,4> Epsd(const xt::xtensor<double,4> &a_Eps);

// compute maximum, avoiding allocation
inline double sigm_max(const xt::xtensor<double,4> &a_Sig);
inline double epsm_max(const xt::xtensor<double,4> &a_Eps);
inline double sigd_max(const xt::xtensor<double,4> &a_Sig);
inline double epsd_max(const xt::xtensor<double,4> &a_Eps);

// ---------------------------- material point - elastic (Elastic.cpp) -----------------------------

class Elastic
{
private:

  // parameters
  double m_kappa; // bulk  modulus
  double m_mu;    // shear modulus

public:

  // constructor
  Elastic() = default;
  Elastic(double kappa, double mu);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double kappa() const;
  double mu() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
};

// -------------------------- material point - cusp potential (Cusp.cpp) ---------------------------

class Cusp
{
private:

  // parameters
  double              m_kappa; // bulk  modulus
  double              m_mu;    // shear modulus
  std::vector<double> m_epsy;  // yield strains

public:

  // constructor
  Cusp() = default;
  Cusp(double kappa, double mu, const std::vector<double> &epsy={}, bool init_elastic=true);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double kappa() const;
  double mu() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
};

// ------------------------ material point - smooth potential (Smooth.cpp) -------------------------

class Smooth
{
private:

  // parameters
  double              m_kappa; // bulk  modulus
  double              m_mu;    // shear modulus
  std::vector<double> m_epsy;  // yield strains

public:

  // constructor
  Smooth() = default;
  Smooth(double kappa, double mu, const std::vector<double> &epsy={}, bool init_elastic=true);

  // stress
  T2s Sig(const T2s &Eps) const;

  // parameters
  double kappa() const;
  double mu() const;

  // energy
  double energy(const T2s &Eps) const;

  // equivalent deviatoric strain
  double epsd(const T2s &Eps) const;

  // index of the current yield strain
  size_t find(const T2s &Eps) const;
  size_t find(double epsd) const;

  // certain yield strain
  double epsy(size_t idx) const;

  // equivalent plastic strain
  double epsp(const T2s &Eps) const;
  double epsp(double epsd) const;
};

// ----------------------- enumerator to switch between material definitions -----------------------

struct Type {
  enum Value {
    Unset,
    Elastic,
    Cusp,
    Smooth,
    PlanarCusp,
    PlanarSmooth,
  };
};

// ------------------- matrix of material points of different types (Matrix.cpp) -------------------

class Matrix
{
private:

  // material vectors
  std::vector<Elastic> m_Elastic;
  std::vector<Cusp>    m_Cusp;
  std::vector<Smooth>  m_Smooth;

  // identifiers for each matrix entry
  xt::xtensor<size_t,2> m_type;  // type (e.g. "Type::Elastic")
  xt::xtensor<size_t,2> m_index; // index from the relevant material vector (e.g. "m_Elastic")

public:

  // constructor
  Matrix() = default;
  Matrix(const std::vector<size_t> &shape);

  // return shape
  std::vector<size_t> shape() const;
  size_t shape(size_t i) const;

  // return type
  xt::xtensor<size_t,2> type() const;

  // return plastic yes/no
  xt::xtensor<size_t,2> isPlastic() const;

  // parameters
  xt::xtensor<double,2> kappa() const;
  xt::xtensor<double,2> mu() const;

  // check that a type has been set everywhere
  void check() const;

  // set material definition for a batch of points
  // -
  void setElastic(const xt::xtensor<size_t,2> &I,
    double kappa, double mu);
  // -
  void setCusp(const xt::xtensor<size_t,2> &I,
    double kappa, double mu, const std::vector<double> &epsy, bool init_elastic=true);
  // -
  void setSmooth(const xt::xtensor<size_t,2> &I,
    double kappa, double mu, const std::vector<double> &epsy, bool init_elastic=true);

  // set material definition for a batch of points
  // -
  void setElastic(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &kappa, const xt::xtensor<double,1> &mu);
  // -
  void setCusp(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &kappa, const xt::xtensor<double,1> &mu,
    const xt::xtensor<double,2> &epsy, bool init_elastic=true);
  // -
  void setSmooth(const xt::xtensor<size_t,2> &I, const xt::xtensor<size_t,2> &idx,
    const xt::xtensor<double,1> &kappa, const xt::xtensor<double,1> &mu,
    const xt::xtensor<double,2> &epsy, bool init_elastic=true);

  // compute (no allocation)
  void Sig   (const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,4> &a_Sig   ) const;
  void energy(const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_energy) const;
  void find  (const xt::xtensor<double,4> &a_Eps, xt::xtensor<size_t,2> &a_find  ) const;
  void epsy  (const xt::xtensor<size_t,2> &a_idx, xt::xtensor<double,2> &a_epsy  ) const;
  void epsp  (const xt::xtensor<double,4> &a_Eps, xt::xtensor<double,2> &a_epsp  ) const;

  // compute (return allocated result)
  xt::xtensor<double,4> Sig   (const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<double,2> energy(const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<size_t,2> find  (const xt::xtensor<double,4> &a_Eps) const;
  xt::xtensor<double,2> epsy  (const xt::xtensor<size_t,2> &a_idx) const;
  xt::xtensor<double,2> epsp  (const xt::xtensor<double,4> &a_Eps) const;

};

// -------------------------------------------------------------------------------------------------

}

// ---------------------------------------- include scripts ----------------------------------------

#include "Cartesian3d.hpp"
#include "Elastic.hpp"
#include "Cusp.hpp"
#include "Smooth.hpp"
#include "Matrix.hpp"

// -------------------------------------------------------------------------------------------------

#endif
