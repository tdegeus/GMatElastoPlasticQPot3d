
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include "../include/GMatElastoPlasticQPot3d/GMatElastoPlasticQPot3d.h"

namespace GM = GMatElastoPlasticQPot3d::Cartesian3d;

// =================================================================================================

TEST_CASE("GMatElastoPlasticQPot3d::Cartesian3d", "Cartesian3d.h")
{

// =================================================================================================

SECTION( "Elastic" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Elastic mat(kappa,mu);

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps *= 0.;
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm  );
  EQ( Sig(1,1), 3. * kappa * epsm  );
  EQ( Sig(2,2), 3. * kappa * epsm  );
  EQ( Sig(0,1), 2. * mu    * gamma );
  EQ( Sig(1,0), 2. * mu    * gamma );
  EQ( Sig(0,2), 0.                 );
  EQ( Sig(2,0), 0.                 );
  EQ( Sig(1,2), 0.                 );
  EQ( Sig(2,1), 0.                 );
  // - plastic strain
  EQ( mat.epsp(Eps), 0 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 0 );
}

// =================================================================================================

SECTION( "Cusp" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Cusp mat(kappa,mu,{0.01, 0.03, 0.10});

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps *= 0.;
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm );
  EQ( Sig(1,1), 3. * kappa * epsm );
  EQ( Sig(2,2), 3. * kappa * epsm );
  EQ( Sig(0,1), 0.                );
  EQ( Sig(1,0), 0.                );
  EQ( Sig(0,2), 0.                );
  EQ( Sig(2,0), 0.                );
  EQ( Sig(1,2), 0.                );
  EQ( Sig(2,1), 0.                );
  // - plastic strain
  EQ( mat.epsp(Eps), 0.02 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 1 );
}

// =================================================================================================

SECTION( "Smooth" )
{
  // material model
  // - parameters
  double kappa = 12.3;
  double mu    = 45.6;
  // - model
  GM::Smooth mat(kappa,mu,{0.01, 0.03, 0.10});

  // allocate tensors
  GM::T2s Eps, Sig;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps *= 0.;
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm );
  EQ( Sig(1,1), 3. * kappa * epsm );
  EQ( Sig(2,2), 3. * kappa * epsm );
  EQ( Sig(0,1), 0.                );
  EQ( Sig(1,0), 0.                );
  EQ( Sig(0,2), 0.                );
  EQ( Sig(2,0), 0.                );
  EQ( Sig(1,2), 0.                );
  EQ( Sig(2,1), 0.                );
  // - plastic strain
  EQ( mat.epsp(Eps), 0.02 );
  // - yield strain index
  REQUIRE( mat.find(Eps) == 1 );
}

// =================================================================================================

SECTION( "Matrix" )
{
  // parameters
  double kappa = 12.3;
  double mu    = 45.6;
  size_t nelem = 3;
  size_t nip   = 2;

  // allocate matrix
  GM::Matrix mat({3,2});

  // row 0: elastic
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>({nelem,nip});

    for ( size_t q = 0 ; q < nip ; ++q ) I(0,q) = 1;

    mat.setElastic(I,kappa,mu);
  }

  // row 1: cups
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>({nelem,nip});

    xt::xtensor<double,1> epsy = {0.01, 0.03, 0.10};

    for ( size_t q = 0 ; q < nip ; ++q ) I(1,q) = 1;

    mat.setCusp(I,kappa,mu,epsy);
  }

  // row 2: smooth
  {
    xt::xtensor<size_t,2> I = xt::zeros<size_t>({nelem,nip});

    xt::xtensor<double,1> epsy = {0.01, 0.03, 0.10};

    for ( size_t q = 0 ; q < nip ; ++q ) I(2,q) = 1;

    mat.setCusp(I,kappa,mu,epsy);
  }

  // allocate tensors
  GM::T2s Eps;

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps *= 0.;
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = Eps(1,0) = gamma;
  // - strain/stress matrices
  xt::xtensor<double,4> eps = xt::empty<double>({3,2,3,3});
  xt::xtensor<double,4> sig;
  xt::xtensor<double,2> epsp;
  // - set strain
  for ( size_t e = 0 ; e < nelem ; ++e ) {
    for ( size_t q = 0 ; q < nip ; ++q ) {
      auto eps_i = xt::view(eps, e, q);
      eps_i = Eps;
    }
  }
  // - stress & plastic strain
  sig  = mat.Sig (eps);
  epsp = mat.epsp(eps);
  // - analytical solution
  EQ( sig(0,0,0,0), 3. * kappa * epsm ); EQ( sig(0,1,0,0), 3. * kappa * epsm );
  EQ( sig(0,0,1,1), 3. * kappa * epsm ); EQ( sig(0,1,1,1), 3. * kappa * epsm );
  EQ( sig(0,0,2,2), 3. * kappa * epsm ); EQ( sig(0,1,2,2), 3. * kappa * epsm );
  EQ( sig(0,0,0,1), 2. * mu    * gamma); EQ( sig(0,1,0,1), 2. * mu    * gamma);
  EQ( sig(0,0,1,0), 2. * mu    * gamma); EQ( sig(0,1,1,0), 2. * mu    * gamma);
  EQ( sig(0,0,0,2), 0.                ); EQ( sig(0,1,0,2), 0.                );
  EQ( sig(0,0,2,0), 0.                ); EQ( sig(0,1,2,0), 0.                );
  EQ( sig(0,0,1,2), 0.                ); EQ( sig(0,1,1,2), 0.                );
  EQ( sig(0,0,2,1), 0.                ); EQ( sig(0,1,2,1), 0.                );
  // -
  EQ( sig(1,0,0,0), 3. * kappa * epsm ); EQ( sig(1,1,0,0), 3. * kappa * epsm );
  EQ( sig(1,0,1,1), 3. * kappa * epsm ); EQ( sig(1,1,1,1), 3. * kappa * epsm );
  EQ( sig(1,0,2,2), 3. * kappa * epsm ); EQ( sig(1,1,2,2), 3. * kappa * epsm );
  EQ( sig(1,0,0,1), 0.                ); EQ( sig(1,1,0,1), 0.                );
  EQ( sig(1,0,1,0), 0.                ); EQ( sig(1,1,1,0), 0.                );
  EQ( sig(1,0,0,2), 0.                ); EQ( sig(1,1,0,2), 0.                );
  EQ( sig(1,0,2,0), 0.                ); EQ( sig(1,1,2,0), 0.                );
  EQ( sig(1,0,1,2), 0.                ); EQ( sig(1,1,1,2), 0.                );
  EQ( sig(1,0,2,1), 0.                ); EQ( sig(1,1,2,1), 0.                );
  // -
  EQ( sig(2,0,0,0), 3. * kappa * epsm ); EQ( sig(2,1,0,0), 3. * kappa * epsm );
  EQ( sig(2,0,1,1), 3. * kappa * epsm ); EQ( sig(2,1,1,1), 3. * kappa * epsm );
  EQ( sig(2,0,2,2), 3. * kappa * epsm ); EQ( sig(2,1,2,2), 3. * kappa * epsm );
  EQ( sig(2,0,0,1), 0.                ); EQ( sig(2,1,0,1), 0.                );
  EQ( sig(2,0,1,0), 0.                ); EQ( sig(2,1,1,0), 0.                );
  EQ( sig(2,0,0,2), 0.                ); EQ( sig(2,1,0,2), 0.                );
  EQ( sig(2,0,2,0), 0.                ); EQ( sig(2,1,2,0), 0.                );
  EQ( sig(2,0,1,2), 0.                ); EQ( sig(2,1,1,2), 0.                );
  EQ( sig(2,0,2,1), 0.                ); EQ( sig(2,1,2,1), 0.                );
  // - plastic strain
  EQ( epsp(0,0), 0    ); EQ( epsp(0,1), 0    );
  EQ( epsp(1,0), gamma); EQ( epsp(1,1), gamma);
  EQ( epsp(2,0), gamma); EQ( epsp(2,1), gamma);
}

// =================================================================================================

}
