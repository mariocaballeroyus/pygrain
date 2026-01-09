"""Export functions for meshes and particle data."""

from __future__ import annotations

import os
from pathlib import Path

import gmsh
import numpy as np
from numpy.typing import NDArray

from ..packing.packing import Packing
from ..packing.particle import Sphere, Spheroid, Cylinder
from ..mesh.mesh import PackingMesh


def to_stl(path: Path | str) -> None:
    """Write the active Gmsh model to a file.
    
    Args:
        path: Output file path (e.g., .msh, .stl, etc.).
    """
    if not gmsh.isInitialized():
        raise RuntimeError(
            "Gmsh is not initialized."
            "Use `Mesh` as a context manager or call `gmsh.initialize()`."
        )
    
    path = Path(path) if isinstance(path, str) else path

    if path.suffix.lower() not in {".msh", ".stl", ".step"}:
        raise ValueError(
            "Unsupported file format. Use: .msh, .stl, .step."
        )
    
    gmsh.write(str(path))
    print(f"\nExported mesh to {path}.")


def to_csv(packing: Packing, path: Path | str, periodic: bool = True) -> None:
    """Export particle data to a single merged CSV file with header sections.
    
    Args:
        packing: Packing object containing particles
        path: Output file path (must end with .csv)
        periodic: Include periodic images if True
    """
    path = Path(path) if isinstance(path, str) else path
    
    if path.suffix.lower() != ".csv":
        raise ValueError("File path must end with .csv extension")
    
    # Create parent directory if it doesn't exist
    if not path.parent.exists():
        os.makedirs(str(path.parent))
    
    data = packing.data_array(periodic)
    indices = np.unique(data[:, 0].astype(int))

    # Storage for each type
    spheres_list = []
    spheroids_list = []
    cylinders_list = []

    for idx in indices:
        if idx >= len(packing.particles):
            continue

        particle_type = packing.particles[idx]
        rows = (data[:, 0].astype(int) == idx)
        subset = data[rows]

        if subset.shape[0] == 0:
            continue

        if isinstance(particle_type, Sphere):
            r = particle_type.radius
            geo = np.full((subset.shape[0], 1), r)
            spheres_list.append(np.hstack([subset[:, 1:4], geo]))

        elif isinstance(particle_type, Spheroid):
            a = particle_type.semi_minor_axis
            c = particle_type.semi_major_axis
            geo = np.full((subset.shape[0], 2), [a, c])
            spheroids_list.append(np.hstack([subset[:, 1:8], geo]))

        elif isinstance(particle_type, Cylinder):
            l = particle_type.length
            d = particle_type.diameter
            geo = np.full((subset.shape[0], 2), [l, d])
            cylinders_list.append(np.hstack([subset[:, 1:8], geo]))

    # Write to single merged CSV file
    total_count = 0
    
    with open(path, 'w') as f:
        # Write spheres section
        if spheres_list:
            f.write("#x,y,z,r\n")
            sphere_data = np.vstack(spheres_list)
            np.savetxt(f, sphere_data, delimiter=",", fmt='%.8f')
            total_count += len(sphere_data)
        
        # Write cylinders section
        if cylinders_list:
            f.write("#x,y,z,ax,ay,az,angle,l,d\n")
            cylinder_data = np.vstack(cylinders_list)
            np.savetxt(f, cylinder_data, delimiter=",", fmt='%.8f')
            total_count += len(cylinder_data)
        
        # Write spheroids section
        if spheroids_list:
            f.write("#x,y,z,ax,ay,az,angle,a,c\n")
            spheroid_data = np.vstack(spheroids_list)
            np.savetxt(f, spheroid_data, delimiter=",", fmt='%.8f')
            total_count += len(spheroid_data)
    
    print(f"Exported {total_count} particles to {path}")
    if periodic:
        print("  * Periodic images included.")