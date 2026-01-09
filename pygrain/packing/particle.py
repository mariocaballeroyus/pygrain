"""Particle type definitions for packing domain."""

from __future__ import annotations

import math
from typing import Protocol, runtime_checkable
from dataclasses import dataclass


@runtime_checkable
class Particle(Protocol):
    """Particle type definition protocol."""
    num: int

    @property
    def volume(self) -> float:
        """Volume of a single particle."""
        ...

    @property
    def acc_volume(self) -> float:
        """Accumulated volume of all particles of this type."""
        ...

    def __eq__(self, other: object) -> bool:
        """Check if two particles are equal (excluding num)."""
        ...


@dataclass(slots=False)
class Sphere(Particle):
    """Sphere particle type definition."""
    num: int

    # Geometric parameters
    radius: float

    # Surface roughness parameters
    corr_length: float = 0.0
    sq_roughness: float = 0.0

    @property
    def volume(self) -> float:
        return (4 / 3) * math.pi * self.radius ** 3
    
    @property
    def acc_volume(self) -> float: 
        return self.num * self.volume
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Sphere):
            return False
        return (abs(self.radius - other.radius) < 1e-10 and
                abs(self.corr_length - other.corr_length) < 1e-10 and
                abs(self.sq_roughness - other.sq_roughness) < 1e-10)


@dataclass(slots=False)
class Spheroid(Particle):
    """Spheroid particle type definition."""
    num: int

    # Geometric parameters
    minor_axis: float
    aspect_ratio: float

    # Surface roughness parameters
    corr_length: float = 0.0
    sq_roughness: float = 0.0

    @property
    def volume(self) -> float:
        return (1 / 3) * math.pi * self.minor_axis**3 * self.aspect_ratio
    
    @property
    def acc_volume(self) -> float: 
        return self.num * self.volume
    
    @property
    def major_axis(self) -> float:
        return self.minor_axis * self.aspect_ratio
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Spheroid):
            return False
        return (abs(self.minor_axis - other.minor_axis) < 1e-10 and 
                abs(self.aspect_ratio - other.aspect_ratio) < 1e-10 and
                abs(self.corr_length - other.corr_length) < 1e-10 and
                abs(self.sq_roughness - other.sq_roughness) < 1e-10)


@dataclass(slots=False)
class Cylinder(Particle):
    """Cylinder particle type definition."""
    num: int

    # Geometric parameters
    diameter: float
    aspect_ratio: float

    # NOTE: 
    # Cylinders do not support surface roughness

    @property
    def volume(self) -> float:
        return (1 / 4) * math.pi * self.diameter**3 * self.aspect_ratio
    
    @property
    def acc_volume(self) -> float: 
        return self.num * self.volume
    
    @property
    def length(self) -> float:
        return self.diameter * self.aspect_ratio
    
    @property
    def radius(self) -> float:
        return self.diameter / 2
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Cylinder):
            return False
        return (abs(self.diameter - other.diameter) < 1e-10 and 
                abs(self.aspect_ratio - other.aspect_ratio) < 1e-10)