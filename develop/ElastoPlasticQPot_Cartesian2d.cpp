
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include "../src/ElastoPlasticQPot3d/ElastoPlasticQPot3d.h"
#include <ElastoPlasticQPot/ElastoPlasticQPot.h>

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
  GM::T2s Eps_GM = GM::T2s::Random();
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
  GM::T2s Eps_GM = GM::T2s::Random();
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
  GM::T2s Eps_GM = GM::T2s::Random();
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
    GM::ArrS I = GM::ArrS::Zero(mat_GM.shape());

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(0,k) = 1;

    mat_GM.setElastic(I,kappa,mu);
  }
  // row 0: elastic
  {
    RF::ArrS I = RF::ArrS::Zero(mat_RF.shape());

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(0,k) = 1;

    mat_RF.setElastic(I,3.*kappa,2.*mu);
  }

  // row 1: cups
  {
    GM::ArrS I = GM::ArrS::Zero(mat_GM.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(1,k) = 1;

    mat_GM.setCusp(I,kappa,mu,epsy);
  }
  // row 1: cups
  {
    RF::ArrS I = RF::ArrS::Zero(mat_RF.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(1,k) = 1;

    mat_RF.setCusp(I,3.*kappa,2.*mu,epsy);
  }

  // row 2: smooth
  {
    GM::ArrS I = GM::ArrS::Zero(mat_GM.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_GM.shape(1) ; ++k ) I(2,k) = 1;

    mat_GM.setCusp(I,kappa,mu,epsy);
  }
  // row 2: smooth
  {
    RF::ArrS I = RF::ArrS::Zero(mat_RF.shape());

    std::vector<double> epsy = {0.01, 0.2, 2.0};

    for ( size_t k = 0 ; k < mat_RF.shape(1) ; ++k ) I(2,k) = 1;

    mat_RF.setCusp(I,3.*kappa,2.*mu,epsy);
  }

  // initialize strain
  // - allocate
  GM::ArrD eps_GM = GM::ArrD::Zero({3,2,GM::T2s::Size()});
  RF::ArrD eps_RF = RF::ArrD::Zero({3,2,RF::T2s::Size()});
  // - fill
  for ( size_t e = 0 ; e < 3 ; ++e ) {
    for ( size_t k = 0 ; k < 2 ; ++k ) {
      // -- random strain
      GM::T2s tmp = GM::T2s::Random();
      // -- store set epsxy
      eps_GM(e,k,1) = tmp(0,0);
      eps_RF(e,k,1) = tmp(0,0);
      // -- store set epsxx
      eps_GM(e,k,0) = tmp(1,1);
      eps_RF(e,k,0) = tmp(1,1);
      // -- store set epsxx
      eps_GM(e,k,3) = -tmp(1,1);
      eps_RF(e,k,2) = -tmp(1,1);
    }
  }

  // - stress & plastic strain
  GM::ArrD sig_GM  = mat_GM.Sig (eps_GM);
  RF::ArrD sig_RF  = mat_RF.Sig (eps_RF);
  GM::ArrD epsp_GM = mat_GM.epsp(eps_GM);
  RF::ArrD epsp_RF = mat_RF.epsp(eps_RF);

  // - check
  for ( size_t e = 0 ; e < 3 ; ++e ) {
    for ( size_t k = 0 ; k < 2 ; ++k ) {
      EQ( sig_GM (e,k,0), sig_RF (e,k,0) );
      EQ( sig_GM (e,k,1), sig_RF (e,k,1) );
      EQ( sig_GM (e,k,3), sig_RF (e,k,2) );
      EQ( epsp_GM(e,k)  , epsp_RF(e,k)   );
      EQ( sig_GM (e,k,2), 0.             );
      EQ( sig_GM (e,k,4), 0.             );
      EQ( sig_GM (e,k,5), 0.             );
    }
  }
}

// =================================================================================================

}
