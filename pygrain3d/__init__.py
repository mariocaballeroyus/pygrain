"""pygrain3d: A Python library for 3D grain packing and mesh generation."""

__version__ = "0.1.0"

from .packing import (
    Packing, 
    create_packing,
    Particle,
    Sphere,
    Spheroid,
    Cylinder
)
from .mesh import PackingMesh
from .io import to_csv, to_stl

__all__ = [
    "Packing",
    "create_packing",
    "Particle",
    "Sphere",
    "Spheroid",
    "Cylinder",
    "PackingMesh",
    "to_csv",
    "to_stl"
]