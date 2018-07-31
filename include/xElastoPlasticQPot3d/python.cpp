/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <pybind11/pybind11.h>
#include <pyxtensor/pyxtensor.hpp>
// #include <xtensor-python/pytensor.hpp>

#include "ElastoPlasticQPot3d.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = ElastoPlasticQPot3d;

// abbreviate types(s)
typedef M::T2s T2s;

// ====================================== ElastoPlasticQPot3d ======================================

PYBIND11_MODULE(ElastoPlasticQPot3d, m) {

m.doc() = "Elasto-plastic material models";

// -------------------------------------------------------------------------------------------------

m.def("epsm", py::overload_cast<const T2s & >(&M::epsm), "Hydrostatic strain" , py::arg("Eps"));
m.def("epsd", py::overload_cast<const T2s & >(&M::epsd), "Eq. strain deviator", py::arg("Eps"));
m.def("Epsd", py::overload_cast<const T2s & >(&M::Epsd), "Strain deviator"    , py::arg("Eps"));

m.def("sigm", py::overload_cast<const T2s & >(&M::sigm), "Hydrostatic stress" , py::arg("Sig"));
m.def("sigd", py::overload_cast<const T2s & >(&M::sigd), "Eq. stress deviator", py::arg("Sig"));
m.def("Sigd", py::overload_cast<const T2s & >(&M::Sigd), "Stress deviator"    , py::arg("Sig"));

m.def("epsm", py::overload_cast<const xt::xtensor<double,4> &>(&M::epsm), "Hydrostatic strain" , py::arg("a_Eps"));
m.def("epsd", py::overload_cast<const xt::xtensor<double,4> &>(&M::epsd), "Eq. strain deviator", py::arg("a_Eps"));
m.def("Epsd", py::overload_cast<const xt::xtensor<double,4> &>(&M::Epsd), "Strain deviator"    , py::arg("a_Eps"));

m.def("sigm", py::overload_cast<const xt::xtensor<double,4> &>(&M::sigm), "Hydrostatic stress" , py::arg("a_Sig"));
m.def("sigd", py::overload_cast<const xt::xtensor<double,4> &>(&M::sigd), "Eq. stress deviator", py::arg("a_Sig"));
m.def("Sigd", py::overload_cast<const xt::xtensor<double,4> &>(&M::Sigd), "Stress deviator"    , py::arg("a_Sig"));

// -------------------------------------------------------------------------------------------------

py::class_<M::Elastic>(m, "Elastic")
  // constructor
  .def(
    py::init<double,double>(),
    "Elastic material",
    py::arg("kappa"),
    py::arg("mu")
  )
  // methods
  .def("Sig"   , &M::Elastic::Sig   , py::arg("Eps"))
  .def("energy", &M::Elastic::energy, py::arg("Eps"))
  .def("epsy"  , &M::Elastic::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&M::Elastic::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&M::Elastic::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&M::Elastic::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&M::Elastic::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const M::Elastic &){
    return "<ElastoPlasticQPot3d.Elastic>"; });

// -------------------------------------------------------------------------------------------------

py::class_<M::Cusp>(m, "Cusp")
  // constructor
  .def(
    py::init<double,double,const std::vector<double>&, bool>(),
    "Cusp material",
    py::arg("kappa"),
    py::arg("mu"),
    py::arg("epsy"),
    py::arg("init_elastic")=true
  )
  // methods
  .def("Sig"   , &M::Cusp::Sig   , py::arg("Eps"))
  .def("energy", &M::Cusp::energy, py::arg("Eps"))
  .def("epsy"  , &M::Cusp::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&M::Cusp::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&M::Cusp::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&M::Cusp::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&M::Cusp::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const M::Cusp &){
    return "<ElastoPlasticQPot3d.Cusp>"; });

// -------------------------------------------------------------------------------------------------

py::class_<M::Smooth>(m, "Smooth")
  // constructor
  .def(
    py::init<double,double,const std::vector<double>&, bool>(),
    "Smooth material",
    py::arg("kappa"),
    py::arg("mu"),
    py::arg("epsy"),
    py::arg("init_elastic")=true
  )
  // methods
  .def("Sig"   , &M::Smooth::Sig   , py::arg("Eps"))
  .def("energy", &M::Smooth::energy, py::arg("Eps"))
  .def("epsy"  , &M::Smooth::epsy  , py::arg("idx"))
  .def("epsp"  , py::overload_cast<const T2s &>(&M::Smooth::epsp, py::const_), py::arg("Eps" ))
  .def("epsp"  , py::overload_cast<double     >(&M::Smooth::epsp, py::const_), py::arg("epsd"))
  .def("find"  , py::overload_cast<const T2s &>(&M::Smooth::find, py::const_), py::arg("Eps" ))
  .def("find"  , py::overload_cast<double     >(&M::Smooth::find, py::const_), py::arg("epsd"))
  // print to screen
  .def("__repr__", [](const M::Smooth &){
    return "<ElastoPlasticQPot3d.Smooth>"; });

// -------------------------------------------------------------------------------------------------

py::module smm = m.def_submodule("Type", "Type enumerator");

py::enum_<M::Type::Value>(smm, "Type")
    .value("Unset"       , M::Type::Unset)
    .value("Elastic"     , M::Type::Elastic)
    .value("Cusp"        , M::Type::Cusp)
    .value("Smooth"      , M::Type::Smooth)
    .value("PlanarCusp"  , M::Type::PlanarCusp)
    .value("PlanarSmooth", M::Type::PlanarSmooth)
    .export_values();

// -------------------------------------------------------------------------------------------------

py::class_<M::Matrix>(m, "Matrix")
  // constructor
  .def(
    py::init<const std::vector<size_t>&>(),
    "Matrix of materials",
    py::arg("shape")
  )
  // methods
  .def("setElastic", py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double                                                            >(&M::Matrix::setElastic),py::arg("I"),               py::arg("kappa"),py::arg("mu"))
  .def("setCusp"   , py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double,                        const std::vector<double>   &, bool>(&M::Matrix::setCusp   ),py::arg("I"),               py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const xt::xtensor<size_t,2> &,                                double,                        double,                        const std::vector<double>   &, bool>(&M::Matrix::setSmooth ),py::arg("I"),               py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setElastic", py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &                                     >(&M::Matrix::setElastic),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"))
  .def("setCusp"   , py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &, const xt::xtensor<double,2> &, bool>(&M::Matrix::setCusp   ),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const xt::xtensor<size_t,2> &, const xt::xtensor<size_t,2> &, const xt::xtensor<double,1> &, const xt::xtensor<double,1> &, const xt::xtensor<double,2> &, bool>(&M::Matrix::setSmooth ),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("shape"     , py::overload_cast<size_t>(&M::Matrix::shape, py::const_))
  .def("shape"     , py::overload_cast<>      (&M::Matrix::shape, py::const_))
  .def("type"      ,                           &M::Matrix::type)
  .def("Sig"       , py::overload_cast<const xt::xtensor<double,4> &>(&M::Matrix::Sig   , py::const_), py::arg("a_Eps"))
  .def("energy"    , py::overload_cast<const xt::xtensor<double,4> &>(&M::Matrix::energy, py::const_), py::arg("a_Eps"))
  .def("find"      , py::overload_cast<const xt::xtensor<double,4> &>(&M::Matrix::find  , py::const_), py::arg("a_Eps"))
  .def("epsy"      , py::overload_cast<const xt::xtensor<size_t,2> &>(&M::Matrix::epsy  , py::const_), py::arg("a_idx"))
  .def("epsp"      , py::overload_cast<const xt::xtensor<double,4> &>(&M::Matrix::epsp  , py::const_), py::arg("a_Eps"))
  // print to screen
  .def("__repr__", [](const M::Matrix &){
    return "<ElastoPlasticQPot3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

