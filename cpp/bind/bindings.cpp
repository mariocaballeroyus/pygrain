#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "packing.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pygrain, m)
{
    py::class_<pygrain::Packing>(m, "Packing")
        .def(py::init<std::array<double, 3>>(), py::arg("lengths"))

        .def("add_sphere_particle", &pygrain::Packing::add_sphere_particle, 
                                    py::arg("radius"))

        .def("add_spheroid_particle", &pygrain::Packing::add_spheroid_particle, 
                                      py::arg("aspect_ratio"), 
                                      py::arg("minor_axis"))

        .def("add_cylinder_particle", &pygrain::Packing::add_cylinder_particle, 
                                      py::arg("aspect_ratio"), 
                                      py::arg("diameter"))

        .def("randomize_particles", &pygrain::Packing::randomize_particles)

        .def("generate", &pygrain::Packing::generate, 
                         py::arg("max_iterations"))

        .def("export_spheres_csv", &pygrain::Packing::export_spheres_csv, 
                                   py::arg("filename"));
}