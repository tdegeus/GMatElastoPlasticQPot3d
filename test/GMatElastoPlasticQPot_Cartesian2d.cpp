
#include <catch2/catch.hpp>

#define EQ(a, b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

#include <GMatElastoPlasticQPot/Cartesian2d.h>
#include <GMatElastoPlasticQPot3d/Cartesian3d.h>

#include <xtensor/xrandom.hpp>

namespace GM = GMatElastoPlasticQPot3d::Cartesian3d;
namespace RF = GMatElastoPlasticQPot::Cartesian2d;

TEST_CASE("GMatElastoPlasticQPot_Cartesian2d", "GMatElastoPlasticQPot_Cartesian2d")
{

double kappa = 12.3;
double mu = 45.6;

GM::Tensor2 Eps_GM = xt::random::rand<double>({3, 3});

// make symmetric
Eps_GM(1, 0) = Eps_GM(0, 1);
Eps_GM(2, 0) = Eps_GM(0, 2);
Eps_GM(2, 1) = Eps_GM(1, 2);

// make plane strain, isochoric
Eps_GM(0, 2) = Eps_GM(2, 0) = Eps_GM(1, 2) = Eps_GM(2, 1) = Eps_GM(2, 2) = 0.;

// make isochoric
Eps_GM(0, 0) = -Eps_GM(1, 1);

// allocate copy
RF::Tensor2 Eps_RF;

// copy
for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
        Eps_RF(i, j) = Eps_GM(i, j);
    }
}

SECTION("Elastic")
{
    GM::Elastic mat_GM(kappa, mu);
    RF::Elastic mat_RF(kappa * 3.0, mu * 2.0);

    GM::Tensor2 Sig_GM = mat_GM.Stress(Eps_GM);
    RF::Tensor2 Sig_RF = mat_RF.Stress(Eps_RF);

    EQ(Sig_GM(0, 0), Sig_RF(0, 0));
    EQ(Sig_GM(1, 1), Sig_RF(1, 1));
    EQ(Sig_GM(0, 1), Sig_RF(0, 1));
    EQ(Sig_GM(1, 0), Sig_RF(1, 0));
    EQ(Sig_GM(0, 2), 0.0);
    EQ(Sig_GM(2, 0), 0.0);
    EQ(Sig_GM(1, 2), 0.0);
    EQ(Sig_GM(2, 1), 0.0);
    EQ(Sig_GM(2, 2), 0.0);

    EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
}

SECTION("Cusp")
{
    GM::Cusp mat_GM(kappa, mu, {0.01, 0.5, 1.0, 2.0});
    RF::Cusp mat_RF(kappa * 3.0, mu * 2.0, {0.01, 0.5, 1.0, 2.0});

    GM::Tensor2 Sig_GM = mat_GM.Stress(Eps_GM);
    RF::Tensor2 Sig_RF = mat_RF.Stress(Eps_RF);

    EQ(Sig_GM(0, 0), Sig_RF(0, 0));
    EQ(Sig_GM(1, 1), Sig_RF(1, 1));
    EQ(Sig_GM(0, 1), Sig_RF(0, 1));
    EQ(Sig_GM(1, 0), Sig_RF(1, 0));
    EQ(Sig_GM(0, 2), 0.);
    EQ(Sig_GM(2, 0), 0.);
    EQ(Sig_GM(1, 2), 0.);
    EQ(Sig_GM(2, 1), 0.);
    EQ(Sig_GM(2, 2), 0.);

    EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
    EQ(mat_GM.epsp(Eps_GM), mat_RF.epsp(Eps_RF));
    EQ(mat_GM.epsy(mat_GM.find(Eps_GM)), mat_RF.epsy(mat_RF.find(Eps_RF)));
    REQUIRE(mat_GM.find(Eps_GM) == mat_RF.find(Eps_RF));
}

SECTION("Smooth")
{
    GM::Smooth mat_GM(kappa, mu, {0.01, 0.5, 1.0, 2.0});
    RF::Smooth mat_RF(kappa * 3.0, mu * 2.0, {0.01, 0.5, 1.0, 2.0});

    GM::Tensor2 Sig_GM = mat_GM.Stress(Eps_GM);
    RF::Tensor2 Sig_RF = mat_RF.Stress(Eps_RF);

    EQ(Sig_GM(0, 0), Sig_RF(0, 0));
    EQ(Sig_GM(1, 1), Sig_RF(1, 1));
    EQ(Sig_GM(0, 1), Sig_RF(0, 1));
    EQ(Sig_GM(1, 0), Sig_RF(1, 0));
    EQ(Sig_GM(0, 2), 0.);
    EQ(Sig_GM(2, 0), 0.);
    EQ(Sig_GM(1, 2), 0.);
    EQ(Sig_GM(2, 1), 0.);
    EQ(Sig_GM(2, 2), 0.);

    EQ(mat_GM.energy(Eps_GM), mat_RF.energy(Eps_RF));
    EQ(mat_GM.epsp(Eps_GM), mat_RF.epsp(Eps_RF));
    EQ(mat_GM.epsy(mat_GM.find(Eps_GM)), mat_RF.epsy(mat_RF.find(Eps_RF)));
    REQUIRE(mat_GM.find(Eps_GM) == mat_RF.find(Eps_RF));
}

SECTION("Matrix")
{
    size_t nelem = 3;
    size_t nip = 2;

    GM::Matrix mat_GM(3, 2);
    RF::Matrix mat_RF(3, 2);

    // row 0: elastic
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_GM.nelem(), mat_GM.nip()});
        xt::view(I, 0, xt::all()) = 1;
        mat_GM.setElastic(I, kappa, mu);
    }
    // row 0: elastic
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_RF.nelem(), mat_RF.nip()});
        xt::view(I, 0, xt::all()) = 1;
        mat_RF.setElastic(I, 3. * kappa, 2. * mu);
    }

    // row 1: cups
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_GM.nelem(), mat_GM.nip()});
        xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
        xt::view(I, 1, xt::all()) = 1;
        mat_GM.setCusp(I, kappa, mu, epsy);
    }
    // row 1: cups
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_RF.nelem(), mat_RF.nip()});
        xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
        xt::view(I, 1, xt::all()) = 1;
        mat_RF.setCusp(I, 3. * kappa, 2. * mu, epsy);
    }

    // row 2: smooth
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_GM.nelem(), mat_GM.nip()});
        xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
        xt::view(I, 2, xt::all()) = 1;
        mat_GM.setCusp(I, kappa, mu, epsy);
    }
    // row 2: smooth
    {
        xt::xtensor<size_t, 2> I = xt::zeros<size_t>({mat_RF.nelem(), mat_RF.nip()});
        xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
        xt::view(I, 2, xt::all()) = 1;
        mat_RF.setCusp(I, 3.0 * kappa, 2.0 * mu, epsy);
    }

    xt::xtensor<double, 4> eps_GM = xt::zeros<double>({nelem, nip, 3ul, 3ul});
    xt::xtensor<double, 4> eps_RF = xt::zeros<double>({nelem, nip, 2ul, 2ul});

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            GM::Tensor2 tmp = xt::random::rand<double>({3, 3});
            eps_GM(e, q, 0, 1) = tmp(0, 1);
            eps_GM(e, q, 1, 0) = tmp(0, 1);
            eps_RF(e, q, 0, 1) = tmp(0, 1);
            eps_RF(e, q, 1, 0) = tmp(0, 1);
            eps_GM(e, q, 0, 0) = tmp(1, 1);
            eps_RF(e, q, 0, 0) = tmp(1, 1);
            eps_GM(e, q, 1, 1) = -tmp(1, 1);
            eps_RF(e, q, 1, 1) = -tmp(1, 1);
        }
    }

    xt::xtensor<double, 4> sig_GM = mat_GM.Stress(eps_GM);
    xt::xtensor<double, 4> sig_RF = mat_RF.Stress(eps_RF);
    xt::xtensor<double, 2> epsp_GM = mat_GM.Epsp(eps_GM);
    xt::xtensor<double, 2> epsp_RF = mat_RF.Epsp(eps_RF);

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            EQ(sig_GM(e, q, 0, 0), sig_RF(e, q, 0, 0));
            EQ(sig_GM(e, q, 0, 1), sig_RF(e, q, 0, 1));
            EQ(sig_GM(e, q, 1, 0), sig_RF(e, q, 1, 0));
            EQ(sig_GM(e, q, 1, 1), sig_RF(e, q, 1, 1));
            EQ(epsp_GM(e, q), epsp_RF(e, q));
            EQ(sig_GM(e, q, 0, 2), 0.0);
            EQ(sig_GM(e, q, 2, 0), 0.0);
            EQ(sig_GM(e, q, 1, 2), 0.0);
            EQ(sig_GM(e, q, 2, 1), 0.0);
            EQ(sig_GM(e, q, 2, 2), 0.0);
        }
    }
}

}
