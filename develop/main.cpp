
#define CATCH_CONFIG_MAIN  // tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT( (a), Catch::WithinAbs((b), 1.e-12) );

#include "../src/ElastoPlasticQPot3d/ElastoPlasticQPot3d.h"

namespace GM = ElastoPlasticQPot3d;

// =================================================================================================

TEST_CASE("ElastoPlasticQPot3d", "ElastoPlasticQPot3d")
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
  GM::T2s Eps = GM::T2s::Zero();
  GM::T2s Sig = GM::T2s::Zero();

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm  );
  EQ( Sig(1,1), 3. * kappa * epsm  );
  EQ( Sig(2,2), 3. * kappa * epsm  );
  EQ( Sig(0,1), 2. * mu    * gamma );
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
  GM::T2s Eps = GM::T2s::Zero();
  GM::T2s Sig = GM::T2s::Zero();

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm );
  EQ( Sig(1,1), 3. * kappa * epsm );
  EQ( Sig(2,2), 3. * kappa * epsm );
  EQ( Sig(0,1), 2. * mu    * 0.   );
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
  GM::T2s Eps = GM::T2s::Zero();
  GM::T2s Sig = GM::T2s::Zero();

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = gamma;
  // - stress
  Sig = mat.Sig(Eps);
  // - analytical solution
  EQ( Sig(0,0), 3. * kappa * epsm );
  EQ( Sig(1,1), 3. * kappa * epsm );
  EQ( Sig(2,2), 3. * kappa * epsm );
  EQ( Sig(0,1), 2. * mu    * 0.   );
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

  // allocate matrix
  GM::Matrix mat({3,2});

  // row 0: elastic
  {
    GM::ArrS I = GM::ArrS::Zero(mat.shape());

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(0,k) = 1;

    mat.setElastic(I,kappa,mu);
  }

  // row 1: cups
  {
    GM::ArrS I = GM::ArrS::Zero(mat.shape());

    std::vector<double> epsy = {0.01, 0.03, 0.10};

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(1,k) = 1;

    mat.setCusp(I,kappa,mu,epsy);
  }

  // row 2: smooth
  {
    GM::ArrS I = GM::ArrS::Zero(mat.shape());

    std::vector<double> epsy = {0.01, 0.03, 0.10};

    for ( size_t k = 0 ; k < mat.shape(1) ; ++k ) I(2,k) = 1;

    mat.setCusp(I,kappa,mu,epsy);
  }

  // allocate tensors
  GM::T2s Eps = GM::T2s::Zero();
  GM::T2s Sig = GM::T2s::Zero();

  // simple shear + volumetric deformation
  // - parameters
  double gamma = 0.02;
  double epsm  = 0.12;
  // - strain
  Eps(0,0) = Eps(1,1) = Eps(2,2) = epsm; Eps(0,1) = gamma;
  // - strain/stress matrices
  GM::ArrD eps({3,2,Eps.size()}), sig, epsp;
  // - set strain
  for ( size_t e = 0 ; e < 3 ; ++e )
    for ( size_t k = 0 ; k < 2 ; ++k )
      std::copy(Eps.begin(), Eps.end(), eps.item(e,k));
  // - stress & plastic strain
  sig  = mat.Sig (eps);
  epsp = mat.epsp(eps);

  // - analytical solution
  EQ( sig(0,0,0), 3. * kappa * epsm ); EQ( sig(0,1,0), 3. * kappa * epsm );
  EQ( sig(0,0,3), 3. * kappa * epsm ); EQ( sig(0,1,3), 3. * kappa * epsm );
  EQ( sig(0,0,5), 3. * kappa * epsm ); EQ( sig(0,1,5), 3. * kappa * epsm );
  EQ( sig(0,0,1), 2. * mu    * gamma); EQ( sig(0,1,1), 2. * mu    * gamma);
  EQ( sig(0,0,2), 2. * mu    * 0.   ); EQ( sig(0,1,2), 2. * mu    * 0.   );
  EQ( sig(0,0,4), 2. * mu    * 0.   ); EQ( sig(0,1,4), 2. * mu    * 0.   );
  // -
  EQ( sig(1,0,0), 3. * kappa * epsm ); EQ( sig(1,1,0), 3. * kappa * epsm );
  EQ( sig(1,0,3), 3. * kappa * epsm ); EQ( sig(1,1,3), 3. * kappa * epsm );
  EQ( sig(1,0,5), 3. * kappa * epsm ); EQ( sig(1,1,5), 3. * kappa * epsm );
  EQ( sig(1,0,1), 2. * mu    * 0.   ); EQ( sig(1,1,1), 2. * mu    * 0.   );
  EQ( sig(1,0,2), 2. * mu    * 0.   ); EQ( sig(1,1,2), 2. * mu    * 0.   );
  EQ( sig(1,0,4), 2. * mu    * 0.   ); EQ( sig(1,1,4), 2. * mu    * 0.   );
  // -
  EQ( sig(2,0,0), 3. * kappa * epsm ); EQ( sig(2,1,0), 3. * kappa * epsm );
  EQ( sig(2,0,3), 3. * kappa * epsm ); EQ( sig(2,1,3), 3. * kappa * epsm );
  EQ( sig(2,0,5), 3. * kappa * epsm ); EQ( sig(2,1,5), 3. * kappa * epsm );
  EQ( sig(2,0,1), 2. * mu    * 0.   ); EQ( sig(2,1,1), 2. * mu    * 0.   );
  EQ( sig(2,0,2), 2. * mu    * 0.   ); EQ( sig(2,1,2), 2. * mu    * 0.   );
  EQ( sig(2,0,4), 2. * mu    * 0.   ); EQ( sig(2,1,4), 2. * mu    * 0.   );
  // - plastic strain
  EQ( epsp(0,0), 0    ); EQ( epsp(0,1), 0    );
  EQ( epsp(1,0), gamma); EQ( epsp(1,1), gamma);
  EQ( epsp(2,0), gamma); EQ( epsp(2,1), gamma);
}

// =================================================================================================

}
