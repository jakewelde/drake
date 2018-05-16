#include "pybind11/eigen.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "drake/bindings/pydrake/pydrake_pybind.h"
#include "drake/bindings/pydrake/util/cpp_template_pybind.h"
#include "drake/bindings/pydrake/systems/systems_pybind.h"
#include "drake/multibody/collision/element.h"
#include "drake/multibody/collision/point_pair.h"
#include "drake/multibody/rigid_body.h"

namespace drake {
namespace pydrake {

PYBIND11_MODULE(collision, m) {
  // NOLINTNEXTLINE(build/namespaces): Emulate placement in namespace.
  using namespace drake::multibody::collision;

  m.doc() = "Drake Collision types.";

  py::module::import("pydrake.multibody.shapes");
  py::module::import("pydrake.multibody.rigid_body");
  py::module::import("pydrake.systems.framework");

   auto bind_common_scalar_types = [m](auto dummy) {
    using T = decltype(dummy);

    DefineTemplateClassWithDefault<PointPair<T>>(
      m, "PointPair", GetPyParam<T>())
      .def_readonly("elementA", &PointPair<T>::elementA, py_reference_internal)
      .def_readonly("elementB", &PointPair<T>::elementB, py_reference_internal)
      .def_readonly("idA", &PointPair<T>::idA)
      .def_readonly("idB", &PointPair<T>::idB)
      .def_readonly("ptA", &PointPair<T>::ptA)
      .def_readonly("ptB", &PointPair<T>::ptB)
      .def_readonly("normal", &PointPair<T>::normal)
      .def_readonly("distance", &PointPair<T>::distance);
  };
  type_visit(bind_common_scalar_types,
             pysystems::CommonScalarPack{});

  py::class_<Element, DrakeShapes::Element>(m, "CollisionElement")
  .def(py::init<const DrakeShapes::Geometry&, const Eigen::Isometry3d&>(),
         py::arg("geometry_in"), py::arg("T_element_to_local"))
  .def("set_body", &Element::set_body)
  .def("get_body", &Element::get_body);
}

}  // namespace pydrake
}  // namespace drake
