from __future__ import annotations

import _pygrain  # type: ignore

from particles import Particle, Sphere, Spheroid, Cylinder


class Packing:
    _cpp_object: _pygrain.Packing

    def __init__(self, lengths: list[float]) -> None:
        self._cpp_object = _pygrain.Packing(lengths)
        self.lengths = lengths
        self.particles: list[Particle] = []

    def add_sphere_particle(self, radius: float) -> None:
        self._cpp_object.add_sphere_particle(radius)
        self.particles.append(Sphere(radius))

    def add_spheroid_particle(self, aspect_ratio: float, minor_axis: float) -> None:
        self._cpp_object.add_spheroid_particle(aspect_ratio, minor_axis)
        self.particles.append(Spheroid(minor_axis, aspect_ratio))

    def add_cylinder_particle(self, aspect_ratio: float, diameter: float) -> None:
        self._cpp_object.add_cylinder_particle(aspect_ratio, diameter)
        self.particles.append(Cylinder(diameter, aspect_ratio))

    def generate(self, max_iterations: int) -> None:
        self._cpp_object.generate(max_iterations)

    @property
    def porosity(self) -> float:
        solid_vol = sum(p.volume for p in self.particles)
        total_vol = self.lengths[0] * self.lengths[1] * self.lengths[2]
        return 1.0 - solid_vol / total_vol


if __name__ == "__main__":
    packing = Packing([100.0, 100.0, 100.0])

    for _ in range(100):
        packing.add_cylinder_particle(2.0, 15.0)

    print(f"Generated packing with porosity = {packing.porosity})")

    print("Generating packing...")
    packing.generate(10000)

    print(f"Packing generated with porosity {packing.porosity}")