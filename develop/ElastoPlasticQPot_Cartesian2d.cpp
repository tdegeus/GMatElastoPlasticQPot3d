
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include "../include/ElastoPlasticQPot3d/ElastoPlasticQPot3d.h"
#include <ElastoPlasticQPot/ElastoPlasticQPot.h>

#include <xtensor/xrandom.hpp>

namespace GM = ElastoPlasticQPot3d;
namespace RF = ElastoPlasticQPot::Cartesian2d;

// =================================================================================================

TEST_CASE("ElastoPlasticQPot_Cartesian2d", "ElastoPlasticQPot_Cartesian2d")
{

// =================================================================================================

SECTION( "Elastic" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Elastic mat_GM(kappa   , mu   );
  RF::Elastic mat_RF(kappa*3., mu*2.);

  // initialize strain
  // - random strain
  GM::T2s Eps_GM = xt::random::rand<double>({GM::ndim, GM::ndim});
  // - make symmetric
  Eps_GM(1,0) = Eps_GM(0,1);
  Eps_GM(2,0) = Eps_GM(0,2);
  Eps_GM(2,1) = Eps_GM(1,2);
  // - make plane strain, isochoric
  Eps_GM(0,2) = Eps_GM(2,0) = Eps_GM(1,2) = Eps_GM(2,1) = Eps_GM(2,2) = 0.;
  // - make isochoric
  Eps_GM(0,0) = -Eps_GM(1,1);
  // - allocate copy
  RF::T2s Eps_RF;
  // copy
  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      Eps_RF(i,j) = Eps_GM(i,j);

  // constitutive response
  GM::T2s Sig_GM = mat_GM.Sig(Eps_GM);
  RF::T2s Sig_RF = mat_RF.Sig(Eps_RF);

  // check
  EQ( Sig_GM(0,0), Sig_RF(0,0) );
  EQ( Sig_GM(1,1), Sig_RF(1,1) );
  EQ( Sig_GM(0,1), Sig_RF(0,1) );
  EQ( Sig_GM(1,0), Sig_RF(1,0) );
  EQ( Sig_GM(0,2), 0.          );
  EQ( Sig_GM(2,0), 0.          );
  EQ( Sig_GM(1,2), 0.          );
  EQ( Sig_GM(2,1), 0.          );
  EQ( Sig_GM(2,2), 0.          );

  // other
  // -
  EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
  EQ(mat_GM.epsd  (Eps_GM), mat_RF.epsd  (Eps_RF));
  EQ(mat_GM.epsp  (Eps_GM), mat_RF.epsp  (Eps_RF));
  // -
  EQ( mat_GM.epsy(mat_GM.find(Eps_GM)) ,  mat_RF.epsy(mat_RF.find(Eps_RF)) );
  // -
  REQUIRE( mat_GM.find(Eps_GM) == mat_RF.find(Eps_RF) );
}

// =================================================================================================

SECTION( "Cusp" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Cusp mat_GM(kappa   , mu   , {0.01, 0.5, 1.0, 2.0});
  RF::Cusp mat_RF(kappa*3., mu*2., {0.01, 0.5, 1.0, 2.0});

  // initialize strain
  // - random strain
  GM::T2s Eps_GM = xt::random::rand<double>({GM::ndim, GM::ndim});
  // - make symmetric
  Eps_GM(1,0) = Eps_GM(0,1);
  Eps_GM(2,0) = Eps_GM(0,2);
  Eps_GM(2,1) = Eps_GM(1,2);
  // - make plane strain, isochoric
  Eps_GM(0,2) = Eps_GM(2,0) = Eps_GM(1,2) = Eps_GM(2,1) = Eps_GM(2,2) = 0.;
  // - make isochoric
  Eps_GM(0,0) = -Eps_GM(1,1);
  // - allocate copy
  RF::T2s Eps_RF;
  // copy
  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      Eps_RF(i,j) = Eps_GM(i,j);

  // constitutive response
  GM::T2s Sig_GM = mat_GM.Sig(Eps_GM);
  RF::T2s Sig_RF = mat_RF.Sig(Eps_RF);

  // check
  EQ( Sig_GM(0,0), Sig_RF(0,0) );
  EQ( Sig_GM(1,1), Sig_RF(1,1) );
  EQ( Sig_GM(0,1), Sig_RF(0,1) );
  EQ( Sig_GM(1,0), Sig_RF(1,0) );
  EQ( Sig_GM(0,2), 0.          );
  EQ( Sig_GM(2,0), 0.          );
  EQ( Sig_GM(1,2), 0.          );
  EQ( Sig_GM(2,1), 0.          );
  EQ( Sig_GM(2,2), 0.          );

  // other
  // -
  EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
  EQ(mat_GM.epsd  (Eps_GM), mat_RF.epsd  (Eps_RF));
  EQ(mat_GM.epsp  (Eps_GM), mat_RF.epsp  (Eps_RF));
  // -
  EQ( mat_GM.epsy(mat_GM.find(Eps_GM)) ,  mat_RF.epsy(mat_RF.find(Eps_RF)) );
  // -
  REQUIRE( mat_GM.find(Eps_GM) == mat_RF.find(Eps_RF) );
}

// =================================================================================================

SECTION( "Smooth" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Smooth mat_GM(kappa   , mu   , {0.01, 0.5, 1.0, 2.0});
  RF::Smooth mat_RF(kappa*3., mu*2., {0.01, 0.5, 1.0, 2.0});

  // initialize strain
  // - random strain
  GM::T2s Eps_GM = xt::random::rand<double>({GM::ndim, GM::ndim});
  // - make symmetric
  Eps_GM(1,0) = Eps_GM(0,1);
  Eps_GM(2,0) = Eps_GM(0,2);
  Eps_GM(2,1) = Eps_GM(1,2);
  // - make plane strain, isochoric
  Eps_GM(0,2) = Eps_GM(2,0) = Eps_GM(1,2) = Eps_GM(2,1) = Eps_GM(2,2) = 0.;
  // - make isochoric
  Eps_GM(0,0) = -Eps_GM(1,1);
  // - allocate copy
  RF::T2s Eps_RF;
  // copy
  for ( size_t i = 0 ; i < 2 ; ++i )
    for ( size_t j = 0 ; j < 2 ; ++j )
      Eps_RF(i,j) = Eps_GM(i,j);

  // constitutive response
  GM::T2s Sig_GM = mat_GM.Sig(Eps_GM);
  RF::T2s Sig_RF = mat_RF.Sig(Eps_RF);

  // check
  EQ( Sig_GM(0,0), Sig_RF(0,0) );
  EQ( Sig_GM(1,1), Sig_RF(1,1) );
  EQ( Sig_GM(0,1), Sig_RF(0,1) );
  EQ( Sig_GM(1,0), Sig_RF(1,0) );
  EQ( Sig_GM(0,2), 0.          );
  EQ( Sig_GM(2,0), 0.          );
  EQ( Sig_GM(1,2), 0.          );
  EQ( Sig_GM(2,1), 0.          );
  EQ( Sig_GM(2,2), 0.          );

  // other
  // -
  EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
  EQ(mat_GM.epsd  (Eps_GM), mat_RF.epsd  (Eps_RF));
  EQ(mat_GM.epsp  (Eps_GM), mat_RF.epsp  (Eps_RF));
  // -
  EQ( mat_GM.epsy(mat_GM.find(Eps_GM)) ,  mat_RF.epsy(mat_RF.find(Eps_RF)) );
  // -
  REQUIRE( mat_GM.find(Eps_GM) == mat_RF.find(Eps_RF) );
}

// =================================================================================================

SECTION( "Matrix" )
{
  // parameters
  double kappa = 12.3;
  double mu    = 45.6;

  // allocate matrix
  GM::Matrix mat_GM({3,2});
  RF::Matrix mat_RF({3,2});

  // row 0: elastic
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_GM.shape());

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(0,k) = 1;

    mat_GM.setElastic(I,kappa,mu);
  }
  // row 0: elastic
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_RF.shape());

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(0,k) = 1;

    mat_RF.setElastic(I,3.*kappa,2.*mu);
  }

  // row 1: cups
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_GM.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(1,k) = 1;

    mat_GM.setCusp(I,kappa,mu,epsy);
  }
  // row 1: cups
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_RF.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(1,k) = 1;

    mat_RF.setCusp(I,3.*kappa,2.*mu,epsy);
  }

  // row 2: smooth
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_GM.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(2,k) = 1;

    mat_GM.setCusp(I,kappa,mu,epsy);
  }
  // row 2: smooth
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>(mat_RF.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(2,k) = 1;

    mat_RF.setCusp(I,3.*kappa,2.*mu,epsy);
  }

  // initialize strain
  // - allocate
  xt::xtensor<double,4> eps_GM = xt::zeros<double>({std::size_t(3),std::size_t(2),GM::ndim,GM::ndim});
  xt::xtensor<double,4> eps_RF = xt::zeros<double>({std::size_t(3),std::size_t(2),RF::ndim,RF::ndim});
  // - fill
  for ( size_t e = 0 ; e < 3 ; ++e ) {
    for ( size_t k = 0 ; k < 2 ; ++k ) {
      // -- random strain
      GM::T2s tmp = xt::random::rand<double>({GM::ndim, GM::ndim});
      // -- store set epsxy
      eps_GM(e,k,0,1) = tmp(0,1);
      eps_GM(e,k,1,0) = tmp(0,1);
      eps_RF(e,k,0,1) = tmp(0,1);
      eps_RF(e,k,1,0) = tmp(0,1);
      // -- store set epsxx
      eps_GM(e,k,0,0) = tmp(1,1);
      eps_RF(e,k,0,0) = tmp(1,1);
      // -- store set epsxx
      eps_GM(e,k,1,1) = -tmp(1,1);
      eps_RF(e,k,1,1) = -tmp(1,1);
    }
  }

  // - stress & plastic strain
  xt::xtensor<double,4> sig_GM  = mat_GM.Sig (eps_GM);
  xt::xtensor<double,4> sig_RF  = mat_RF.Sig (eps_RF);
  xt::xtensor<double,2> epsp_GM = mat_GM.epsp(eps_GM);
  xt::xtensor<double,2> epsp_RF = mat_RF.epsp(eps_RF);

  // - check
  for ( size_t e = 0 ; e < 3 ; ++e ) {
    for ( size_t k = 0 ; k < 2 ; ++k ) {
      EQ( sig_GM (e,k,0,0), sig_RF (e,k,0,0) );
      EQ( sig_GM (e,k,0,1), sig_RF (e,k,0,1) );
      EQ( sig_GM (e,k,1,0), sig_RF (e,k,1,0) );
      EQ( sig_GM (e,k,1,1), sig_RF (e,k,1,1) );
      EQ( epsp_GM(e,k)    , epsp_RF(e,k)     );
      EQ( sig_GM (e,k,0,2), 0.               );
      EQ( sig_GM (e,k,2,0), 0.               );
      EQ( sig_GM (e,k,1,2), 0.               );
      EQ( sig_GM (e,k,2,1), 0.               );
      EQ( sig_GM (e,k,2,2), 0.               );
    }
  }
}

// =================================================================================================

}
