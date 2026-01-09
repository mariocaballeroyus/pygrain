"""Packing domain management with multiple particles."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .. import _pygrain

from .particle import Particle, Sphere, Spheroid, Cylinder


class Packing:
    """Packing domain containing multiple particles."""
    _cpp_object: _pygrain.Packing

    def __init__(self, lengths: list[float]) -> None:
        """Initialize a packing domain.
        
        Args:
            lengths: Domain lengths in x, y, z directions.
        """
        self._cpp_object = _pygrain.Packing(lengths)
        self.lengths = lengths
        self.particles: list[Particle] = []

    def add_sphere_particles(self, radius: float, num: int, corr_length: float = 0.0, sq_roughness: float = 0.0) -> int:
        """Add multiple sphere particles to the packing.
        
        Args:
            radius: Sphere radius.
            num: Number of sphere particles to add.
            corr_length: Correlation length for surface roughness (default: 0.0, no roughness).
            sq_roughness: Squared roughness parameter (variance) (default: 0.0, no roughness).
            
        Returns:
            The geometry index (ID) assigned to these particles.
        """
        sphere = Sphere(num=num, radius=radius, corr_length=corr_length, sq_roughness=sq_roughness)

        for idx, p in enumerate(self.particles):
            if p == sphere:
                p.num += num
                self._cpp_object.add_sphere_particles(radius, num, idx)
                return idx
        
        idx = len(self.particles)
        self.particles.append(sphere)
        self._cpp_object.add_sphere_particles(radius, num, idx)
        return idx

    def add_spheroid_particles(self, aspect_ratio: float, minor_axis: float, num: int, corr_length: float = 0.0, sq_roughness: float = 0.0) -> int:
        """Add multiple spheroid particles to the packing.
        
        Args:
            aspect_ratio: Spheroid aspect ratio (major_axis / minor_axis).
            minor_axis: Minor axis length (full axis).
            num: Number of spheroid particles to add.
            corr_length: Correlation length for surface roughness (default: 0.0, no roughness).
            sq_roughness: Squared roughness parameter (variance) (default: 0.0, no roughness).
            
        Returns:
            The geometry index (ID) assigned to these particles.
        """
        spheroid = Spheroid(num=num, minor_axis=minor_axis, aspect_ratio=aspect_ratio, corr_length=corr_length, sq_roughness=sq_roughness)

        for idx, p in enumerate(self.particles):
            if p == spheroid:
                p.num += num
                self._cpp_object.add_spheroid_particles(aspect_ratio, minor_axis, num, idx)
                return idx
        
        idx = len(self.particles)
        self.particles.append(spheroid)
        self._cpp_object.add_spheroid_particles(aspect_ratio, minor_axis, num, idx)
        return idx

    def add_cylinder_particles(self, aspect_ratio: float, diameter: float, num: int) -> int:
        """Add multiple cylinder particles to the packing.
        
        Args:
            aspect_ratio: Cylinder aspect ratio (height / diameter).
            diameter: Cylinder diameter (full diameter).
            num: Number of cylinder particles to add.
            
        Returns:
            The geometry index (ID) assigned to these particles.
        """
        cylinder = Cylinder(num=num, diameter=diameter, aspect_ratio=aspect_ratio)

        for idx, p in enumerate(self.particles):
            if p == cylinder:
                p.num += num
                self._cpp_object.add_cylinder_particles(aspect_ratio, diameter, num, idx)
                return idx
        
        idx = len(self.particles)
        self.particles.append(cylinder)
        self._cpp_object.add_cylinder_particles(aspect_ratio, diameter, num, idx)
        return idx

    def generate(self, max_iterations: int, log_interval: int = 500) -> None:
        """Generate a packing with non-overlapping particles.
        
        Args:
            max_iterations: Maximum number of iterations to attempt packing.
            log_interval: Interval for logging progress.
        """
        print(f"Starting packing generation with {self.num_particles} particles ...")
        self._cpp_object.generate(max_iterations, log_interval)

    @property
    def num_particles(self) -> int:
        """Get the total number of particles in the packing."""
        return self._cpp_object.num_particles()

    @property
    def volume(self) -> float:
        """Get the total volume of the packing."""
        return self.lengths[0] * self.lengths[1] * self.lengths[2]

    @property
    def porosity(self) -> float:
        """Calculate and return the porosity of the packing."""
        solid_vol = sum(p.acc_volume for p in self.particles)
        return 1.0 - solid_vol / self.volume

    def data_array(self, periodic: bool = False) -> NDArray[np.float64]:
        """Get all particle positions a numpy array.
        
        Args:
            periodic: Whether to export the peridic images of particles.

        Returns:
            Numpy array of shape (num_particles, 8) with columns:
            [idx, x, y, z, axis_x, axis_y, axis_z, angle]
        """
        return self._cpp_object.data_array(periodic)
    

def create_packing(lengths: list[float]) -> Packing:
    """Create a packing domain.
    
    Args:
        lengths: Domain lengths in x, y, z directions.

    Returns:
        Packing object.
    """
    if len(lengths) != 3:
        raise ValueError("Lengths must be a list of three dimensions (x, y, z).")
    
    return Packing(lengths)