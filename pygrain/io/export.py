from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import gmsh

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
    """Export particle data to a CSV file. The file will contain separate 
    chunks for each supported particle type.
    
    Args:
        packing: Packing instance containing particle data.
        path: Output CSV file path.
        periodic: Whether to include periodic images (default: True).
    """
    path = Path(path) if isinstance(path, str) else path
    
    if path.suffix.lower() != ".csv":
        raise ValueError("Output file must have a .csv extension.")
    
    data = packing.data_array(periodic)
    indices = np.unique(data[:, 0].astype(int))

    spheres_data: list[NDArray[np.float64]]  = []
    spheroids_data: list[NDArray[np.float64]] = []
    cylinders_data: list[NDArray[np.float64]] = []

    for idx in indices:
        if idx >= len(packing.particles):
            continue  # skip empty indices

        particle_type = packing.particles[idx]
        rows = (data[:, 0].astype(int) == idx)
        subset = data[rows]

        if subset.shape[0] == 0:
            continue  # skip if no data for this particle

        if isinstance(particle_type, Sphere):
            # Format: [x, y, z, r]
            r = particle_type.radius
            geo = np.full((subset.shape[0], 1), r)  # r
            pos = subset[:, 1:4]                    # x, y, z

            chunk = np.hstack([pos, geo])
            spheres_data.append(chunk)

        elif isinstance(particle_type, Spheroid):
            # Format: [x, y, z, axis_x, axis_y, axis_z, angle, a, c]
            a = particle_type.semi_minor_axis
            c = particle_type.semi_major_axis
            geo = np.full((subset.shape[0], 2), [a, c])  # a, c
            pos = subset[:, 1:8]                         # x, y, z, ax, ay, az, angle

            chunk = np.hstack([pos, geo])
            spheroids_data.append(chunk)

        elif isinstance(particle_type, Cylinder):
            # Format: [x, y, z, axis_x, axis_y, axis_z, angle, L, D]
            l = particle_type.length
            d = particle_type.diameter
            geo = np.full((subset.shape[0], 2), [l, d])  # l, d
            pos = subset[:, 1:8]                         # x, y, z, ax, ay, az, angle

            chunk = np.hstack([pos, geo])
            cylinders_data.append(chunk)

    with open(path, 'w') as f:
        if spheres_data:
            header = "x,y,z,r"
            all_spheres = np.vstack(spheres_data)
            np.savetxt(f, all_spheres, delimiter=",", header=header, fmt='%.8f')

        if spheroids_data:
            header = "x,y,z,axis_x,axis_y,axis_z,angle,a,c"
            all_spheroids = np.vstack(spheroids_data)
            np.savetxt(f, all_spheroids, delimiter=",", header=header, fmt='%.8f')

        if cylinders_data:
            header = "x,y,z,axis_x,axis_y,axis_z,angle,L,D"
            all_cylinders = np.vstack(cylinders_data)
            np.savetxt(f, all_cylinders, delimiter=",", header=header, fmt='%.8f')

    print(f"\nExported particle data to {path}.")
    if periodic:
        print("  * Periodic images included.")