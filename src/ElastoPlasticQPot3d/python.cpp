/* =================================================================================================

(c - GPLv3) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | github.com/tdegeus/GooseFEM

================================================================================================= */

#include <Eigen/Eigen>
#include <cppmat/cppmat.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <cppmat/pybind11.h>

#include "ElastoPlasticQPot3d.h"

// =================================================================================================

// abbreviate name-space
namespace py = pybind11;
namespace M  = ElastoPlasticQPot3d;

// abbreviate types(s)
typedef M::T2s  T2s;
typedef M::ArrD ArrD;
typedef M::ArrS ArrS;
typedef M::MatD MatD;
typedef M::ColD ColD;

// ====================================== ElastoPlasticQPot3d ======================================

PYBIND11_MODULE(ElastoPlasticQPot3d, m) {

m.doc() = "Elasto-plastic material models";

// -------------------------------------------------------------------------------------------------

m.def("epsm", py::overload_cast<const T2s & >(&M::epsm), "Mean strain"        , py::arg("Eps"));
m.def("epsd", py::overload_cast<const T2s & >(&M::epsd), "Eq. strain deviator", py::arg("Eps"));
m.def("Epsd", py::overload_cast<const T2s & >(&M::Epsd), "Strain deviator"    , py::arg("Eps"));

m.def("sigm", py::overload_cast<const T2s & >(&M::sigm), "Mean stress"        , py::arg("Sig"));
m.def("sigd", py::overload_cast<const T2s & >(&M::sigd), "Eq. stress deviator", py::arg("Sig"));
m.def("Sigd", py::overload_cast<const T2s & >(&M::Sigd), "Stress deviator"    , py::arg("Sig"));

m.def("epsm", py::overload_cast<const ArrD &>(&M::epsm), "Mean strain"        , py::arg("a_Eps"));
m.def("epsd", py::overload_cast<const ArrD &>(&M::epsd), "Eq. strain deviator", py::arg("a_Eps"));
m.def("Epsd", py::overload_cast<const ArrD &>(&M::Epsd), "Strain deviator"    , py::arg("a_Eps"));

m.def("sigm", py::overload_cast<const ArrD &>(&M::sigm), "Mean stress"        , py::arg("a_Sig"));
m.def("sigd", py::overload_cast<const ArrD &>(&M::sigd), "Eq. stress deviator", py::arg("a_Sig"));
m.def("Sigd", py::overload_cast<const ArrD &>(&M::Sigd), "Stress deviator"    , py::arg("a_Sig"));

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
  .def("__repr__", [](const M::Elastic &a){
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
  .def("__repr__", [](const M::Cusp &a){
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
  .def("__repr__", [](const M::Smooth &a){
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
  .def("setElastic", py::overload_cast<const ArrS &, double, double                                              >(&M::Matrix::setElastic),py::arg("I"),               py::arg("kappa"),py::arg("mu"))
  .def("setCusp"   , py::overload_cast<const ArrS &, double, double, const std::vector<double> &, bool           >(&M::Matrix::setCusp   ),py::arg("I"),               py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const ArrS &, double, double, const std::vector<double> &, bool           >(&M::Matrix::setSmooth ),py::arg("I"),               py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setElastic", py::overload_cast<const ArrS &, const ArrS &, const ColD &, const ColD &                    >(&M::Matrix::setElastic),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"))
  .def("setCusp"   , py::overload_cast<const ArrS &, const ArrS &, const ColD &, const ColD &, const MatD &, bool>(&M::Matrix::setCusp   ),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("setSmooth" , py::overload_cast<const ArrS &, const ArrS &, const ColD &, const ColD &, const MatD &, bool>(&M::Matrix::setSmooth ),py::arg("I"),py::arg("idx"),py::arg("kappa"),py::arg("mu"),py::arg("epsy"),py::arg("init_elastic")=true)
  .def("shape"     , py::overload_cast<size_t>(&M::Matrix::shape, py::const_))
  .def("shape"     , py::overload_cast<      >(&M::Matrix::shape, py::const_))
  .def("type"      , &M::Matrix::type)
  .def("Sig"       , &M::Matrix::Sig   , py::arg("a_Eps"))
  .def("energy"    , &M::Matrix::energy, py::arg("a_Eps"))
  .def("epsy"      , &M::Matrix::epsy  , py::arg("a_idx"))
  .def("epsp"      , &M::Matrix::epsp  , py::arg("a_Eps"))
  .def("find"      , &M::Matrix::find  , py::arg("a_Eps"))
  // print to screen
  .def("__repr__", [](const M::Matrix &a){
    return "<ElastoPlasticQPot3d.Matrix>"; });

// -------------------------------------------------------------------------------------------------

}

