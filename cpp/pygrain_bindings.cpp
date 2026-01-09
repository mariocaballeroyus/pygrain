#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "packing.hpp"

namespace py = pybind11;

namespace pygrain3d
{

PYBIND11_MODULE(_pygrain3d, m)
{
    py::class_<Packing>(m, "Packing")
        .def(py::init<std::array<double, 3>>(), py::arg("lengths"))

        .def("add_sphere_particles", &Packing::add_sphere_particles, 
                                    py::arg("radius"),
                                    py::arg("num"),
                                    py::arg("id"))

        .def("add_spheroid_particles", &Packing::add_spheroid_particles, 
                                      py::arg("aspect_ratio"), 
                                      py::arg("minor_axis"),
                                      py::arg("num"),
                                      py::arg("id"))

        .def("add_cylinder_particles", &Packing::add_cylinder_particles, 
                                      py::arg("aspect_ratio"), 
                                      py::arg("diameter"),
                                      py::arg("num"),
                                      py::arg("id"))

        .def("randomize_particles", &Packing::randomize_particles)

        .def("generate", &Packing::generate, 
                         py::arg("max_iterations"),
                         py::arg("log_interval"))

        .def("num_particles", &Packing::num_particles)
        .def("clear", &Packing::clear)

        .def("particle_array", [](const Packing& self, bool periodic) {
            auto positions = self.particle_array(periodic);
            std::size_t n = positions.size() / 8;
            return py::array_t<double>({n, std::size_t(8)}, positions.data());
        }, py::arg("periodic"))

        .def("sphere_array", [](const Packing& self) {
            auto spheres = self.sphere_array();
            std::size_t n = spheres.size() / 4;
            return py::array_t<double>({n, std::size_t(4)}, spheres.data());
        });
}

} // namespace pygrain3d