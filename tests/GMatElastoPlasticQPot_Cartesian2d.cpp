#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <xtensor/xrandom.hpp>
#include <GMatElastoPlasticQPot/Cartesian2d.h>
#include <GMatElastoPlasticQPot3d/Cartesian3d.h>

#define ISCLOSE(a, b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1e-12));

namespace GM = GMatElastoPlasticQPot3d::Cartesian3d;
namespace RF = GMatElastoPlasticQPot::Cartesian2d;

template <class T, class S>
S A4_ddot_B2(const T& A, const S& B)
{
    S C = xt::empty<double>({3, 3});
    C.fill(0.0);

    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < 3; k++) {
                for (size_t l = 0; l < 3; l++) {
                    C(i, j) += A(i, j, k, l) * B(l, k);
                }
            }
        }
    }

    return C;
}

TEST_CASE("GMatElastoPlasticQPot_Cartesian2d", "GMatElastoPlasticQPot_Cartesian2d")
{
    SECTION("Elastic")
    {
        double kappa = 12.3;
        double mu = 45.6;

        // symmetric random strain
        xt::xtensor<double, 2> Eps_GM = xt::random::rand<double>({3, 3});
        auto I4s = GM::I4s();
        Eps_GM = A4_ddot_B2(I4s, Eps_GM);

        // make plane strain
        Eps_GM(0, 2) = Eps_GM(2, 0) = Eps_GM(1, 2) = Eps_GM(2, 1) = Eps_GM(2, 2) = 0.0;

        // make isochoric
        Eps_GM(0, 0) = -Eps_GM(1, 1);

        // copy
        xt::xtensor<double, 2> Eps_RF = xt::view(Eps_GM, xt::range(0, 2), xt::range(0, 2));

        GM::Elastic mat_GM(kappa, mu);
        RF::Elastic mat_RF(kappa * 3.0, mu * 2.0);

        mat_GM.setStrain(Eps_GM);
        mat_RF.setStrain(Eps_RF);

        auto Sig_GM = mat_GM.Stress();
        auto Sig_RF = mat_RF.Stress();

        xt::xtensor<double, 2> Sig = xt::zeros<double>({3, 3});
        xt::view(Sig, xt::range(0, 2), xt::range(0, 2)) = Sig_RF;

        REQUIRE(xt::allclose(Sig, Sig_GM));
        REQUIRE(mat_GM.energy() == Approx(mat_RF.energy()));
    }

    SECTION("Cusp")
    {
        double kappa = 12.3;
        double mu = 45.6;

        // symmetric random strain
        xt::xtensor<double, 2> Eps_GM = xt::random::rand<double>({3, 3});
        auto I4s = GM::I4s();
        Eps_GM = A4_ddot_B2(I4s, Eps_GM);

        // make plane strain
        Eps_GM(0, 2) = Eps_GM(2, 0) = Eps_GM(1, 2) = Eps_GM(2, 1) = Eps_GM(2, 2) = 0.0;

        // make isochoric
        Eps_GM(0, 0) = -Eps_GM(1, 1);

        // copy
        xt::xtensor<double, 2> Eps_RF = xt::view(Eps_GM, xt::range(0, 2), xt::range(0, 2));

        GM::Cusp mat_GM(kappa, mu, {0.01, 0.5, 1.0, 2.0});
        RF::Cusp mat_RF(kappa * 3.0, mu * 2.0, xt::xtensor<double , 1>{0.01, 0.5, 1.0, 2.0});

        mat_GM.setStrain(Eps_GM);
        mat_RF.setStrain(Eps_RF);

        auto Sig_GM = mat_GM.Stress();
        auto Sig_RF = mat_RF.Stress();

        xt::xtensor<double, 2> Sig = xt::zeros<double>({3, 3});
        xt::view(Sig, xt::range(0, 2), xt::range(0, 2)) = Sig_RF;

        REQUIRE(xt::allclose(Sig, Sig_GM));
        REQUIRE(mat_GM.energy() == Approx(mat_RF.energy()));
        REQUIRE(mat_GM.epsp() == Approx(mat_RF.epsp()));
        REQUIRE(mat_GM.currentYieldLeft() == Approx(mat_RF.currentYieldLeft()));
        REQUIRE(mat_GM.currentYieldRight() == Approx(mat_RF.currentYieldRight()));
        REQUIRE(static_cast<long>(mat_GM.currentIndex()) == mat_RF.currentIndex());
    }

    SECTION("Smooth")
    {
        double kappa = 12.3;
        double mu = 45.6;

        // symmetric random strain
        xt::xtensor<double, 2> Eps_GM = xt::random::rand<double>({3, 3});
        auto I4s = GM::I4s();
        Eps_GM = A4_ddot_B2(I4s, Eps_GM);

        // make plane strain
        Eps_GM(0, 2) = Eps_GM(2, 0) = Eps_GM(1, 2) = Eps_GM(2, 1) = Eps_GM(2, 2) = 0.0;

        // make isochoric
        Eps_GM(0, 0) = -Eps_GM(1, 1);

        // copy
        xt::xtensor<double, 2> Eps_RF = xt::view(Eps_GM, xt::range(0, 2), xt::range(0, 2));

        GM::Smooth mat_GM(kappa, mu, {0.01, 0.5, 1.0, 2.0});
        RF::Smooth mat_RF(kappa * 3.0, mu * 2.0, xt::xtensor<double , 1>{0.01, 0.5, 1.0, 2.0});

        mat_GM.setStrain(Eps_GM);
        mat_RF.setStrain(Eps_RF);

        auto Sig_GM = mat_GM.Stress();
        auto Sig_RF = mat_RF.Stress();

        xt::xtensor<double, 2> Sig = xt::zeros<double>({3, 3});
        xt::view(Sig, xt::range(0, 2), xt::range(0, 2)) = Sig_RF;

        REQUIRE(xt::allclose(Sig, Sig_GM));
        REQUIRE(mat_GM.energy() == Approx(mat_RF.energy()));
        REQUIRE(mat_GM.epsp() == Approx(mat_RF.epsp()));
        REQUIRE(mat_GM.currentYieldLeft() == Approx(mat_RF.currentYieldLeft()));
        REQUIRE(mat_GM.currentYieldRight() == Approx(mat_RF.currentYieldRight()));
        REQUIRE(static_cast<long>(mat_GM.currentIndex()) == mat_RF.currentIndex());
    }

    SECTION("Array")
    {
        double kappa = 12.3;
        double mu = 45.6;

        size_t nelem = 3;
        size_t nip = 2;

        GM::Array<2> mat_GM({nelem, nip});
        RF::Array<2> mat_RF({nelem, nip});

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, 0, xt::all()) = 1;
            mat_GM.setElastic(I, kappa, mu);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::view(I, 0, xt::all()) = 1;
            mat_RF.setElastic(I, 3. * kappa, 2. * mu);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
            xt::view(I, 1, xt::all()) = 1;
            mat_GM.setCusp(I, kappa, mu, epsy);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
            xt::view(I, 1, xt::all()) = 1;
            mat_RF.setCusp(I, 3. * kappa, 2. * mu, epsy);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
            xt::view(I, 2, xt::all()) = 1;
            mat_GM.setCusp(I, kappa, mu, epsy);
        }

        {
            xt::xtensor<size_t, 2> I = xt::zeros<size_t>({nelem, nip});
            xt::xtensor<double, 1> epsy = {0.01, 0.2, 2.0};
            xt::view(I, 2, xt::all()) = 1;
            mat_RF.setCusp(I, 3.0 * kappa, 2.0 * mu, epsy);
        }

        xt::xtensor<double, 4> eps_GM = xt::zeros<double>({nelem, nip, size_t(3), size_t(3)});
        xt::xtensor<double, 4> eps_RF = xt::zeros<double>({nelem, nip, size_t(2), size_t(2)});

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                xt::xtensor<double, 2> tmp = xt::random::rand<double>({3, 3});
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

        mat_GM.setStrain(eps_GM);
        mat_RF.setStrain(eps_RF);

        auto sig_GM = mat_GM.Stress();
        auto sig_RF = mat_RF.Stress();
        auto epsp_GM = mat_GM.Epsp();
        auto epsp_RF = mat_RF.Epsp();

        for (size_t e = 0; e < nelem; ++e) {
            for (size_t q = 0; q < nip; ++q) {
                ISCLOSE(sig_GM(e, q, 0, 0), sig_RF(e, q, 0, 0));
                ISCLOSE(sig_GM(e, q, 0, 1), sig_RF(e, q, 0, 1));
                ISCLOSE(sig_GM(e, q, 1, 0), sig_RF(e, q, 1, 0));
                ISCLOSE(sig_GM(e, q, 1, 1), sig_RF(e, q, 1, 1));
                ISCLOSE(epsp_GM(e, q), epsp_RF(e, q));
                ISCLOSE(sig_GM(e, q, 0, 2), 0.0);
                ISCLOSE(sig_GM(e, q, 2, 0), 0.0);
                ISCLOSE(sig_GM(e, q, 1, 2), 0.0);
                ISCLOSE(sig_GM(e, q, 2, 1), 0.0);
                ISCLOSE(sig_GM(e, q, 2, 2), 0.0);
            }
        }
    }

}
