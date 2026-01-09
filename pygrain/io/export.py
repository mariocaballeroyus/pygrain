from __future__ import annotations

from pathlib import Path

import numpy as np
from numpy.typing import NDArray

import gmsh

from ..packing.packing import Packing
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
    """Export particle data to a CSV file.
    
    Args:
        packing: Packing instance containing particle data.
        filename: Output CSV file path.
        periodic: Whether to include periodic images (default: True).
    """
    path = Path(path) if isinstance(path, str) else path
    
    if path.suffix.lower() != ".csv":
        raise ValueError("Output file must have a .csv extension.")
    
    data: NDArray[np.float64] = packing.data_array(periodic)
    header = "id,x,y,z,axis_x,axis_y,axis_z,angle"
    format_list = ['%d',                    # id
                   '%.8f', '%.8f', '%.8f',  # x, y, z
                   '%.8f', '%.8f', '%.8f',  # axis_x, axis_y, axis_z
                   '%.8f']                  # angle

    np.savetxt(
        path, data, delimiter=",", header=header, fmt=format_list
    )

    print(f"\nExported {len(data)} particles to {path}.")
    if periodic:
        print("  * Periodic images included.")