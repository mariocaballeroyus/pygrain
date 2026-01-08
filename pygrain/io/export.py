from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

import gmsh

if TYPE_CHECKING:
    from ..mesh import Mesh

from ..packing.packing import Packing


def to_stl(mesh: "Mesh", path: str) -> None:
    """Write a mesh to a file.
    
    Args:
        mesh: The Mesh object to write.
        path: Output file path. Format determined by extension:
              .msh (Gmsh native), .stl, .vtk, .inp (Abaqus), etc.
    
    Example:
        with Mesh("packing") as mesh:
            mesh.generate_from_packing(packing, mesh_size=2.0)
            to_stl(mesh, "output.stl")
    """
    if not gmsh.isInitialized():
        raise RuntimeError(
            "Gmsh is not initialized. "
            "Use `Mesh` as a context manager or call gmsh.initialize()"
        )
    gmsh.write(path)


def to_csv(packing: Packing, filename: str, periodic: bool = True) -> None:
    data = packing.data_array(periodic)
    
    header = "id,x,y,z,axis_x,axis_y,axis_z,angle"
    format_list = ['%d',                    # id
                   '%.8f', '%.8f', '%.8f',  # x, y, z
                   '%.8f', '%.8f', '%.8f',  # axis_x, axis_y, axis_z
                   '%.8f']                  # angle

    np.savetxt(
        filename, data, delimiter=",", header=header, fmt=format_list
    )

    print(f"Exported {len(data)} particles to {filename}.")
    if periodic:
        print("  (including periodic images)")