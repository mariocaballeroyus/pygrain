from __future__ import annotations

import gmsh
import numpy as np
from numpy.typing import NDArray
from scipy.spatial.transform import Rotation

from ..packing.particle import Particle, Sphere, Spheroid, Cylinder


def _create_occ_geometry(particle: Particle) -> int:
    """Create a template geometry at the origin for a particle type.
    
    The template is created at the origin with the reference orientation
    (aligned along the x-axis).
    
    Args:
        particle: The particle type definition.
    
    Returns:
        The gmsh volume tag of the template geometry.
    """
    if isinstance(particle, Sphere):
        return gmsh.model.occ.addSphere(0, 0, 0, particle.radius)
    
    elif isinstance(particle, Spheroid):
        semi_minor = particle.minor_axis / 2.0
        semi_major = particle.major_axis / 2.0
        vol = gmsh.model.occ.addSphere(0, 0, 0, 1)
        gmsh.model.occ.dilate([(3, vol)], 0, 0, 0, semi_major, semi_minor, semi_minor)
        return vol
    
    elif isinstance(particle, Cylinder):
        length = particle.length
        radius = particle.radius
        return gmsh.model.occ.addCylinder(-length / 2, 0, 0, length, 0, 0, radius)
    
    # TODO: 
    # Add support for new particles in the future
    
    else:
        raise TypeError(f"Unknown particle type: {type(particle)}")


def _mesh_occ_geometry(
    mesh_size: float
) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
    """Create OCC geometry, mesh it, and return nodes and elements.
    
    Args:
        mesh_size: Target mesh element size.
    
    Returns:
        Tuple of (nodes, triangles) as numpy arrays.
    """
    gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size * 0.5)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size)
    gmsh.model.mesh.generate(2)

    ntags, coords, _ = gmsh.model.mesh.getNodes()
    nodes = coords.reshape(-1, 3)

    etypes, _, elem_ntags = gmsh.model.mesh.getElements(dim=2)

    if len(etypes) == 0:
        raise ValueError("No 2D elements found in the mesh")

    idx = list(etypes).index(2)  # 2D3N elements
    tri_node_tags = elem_ntags[idx].reshape(-1, 3)

    sort_idx = np.argsort(ntags)
    elements = np.searchsorted(ntags, tri_node_tags, sorter=sort_idx)
    elements = sort_idx[elements]
    
    return nodes, elements


def _transform_nodes(
    nodes: NDArray[np.float64],
    axis: NDArray[np.float64],
    angle: float,
    translation: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Transform nodes by rotation then translation.
    
    Args:
        nodes: Array of node coordinates (N x 3).
        axis: Rotation axis (3,).
        angle: Rotation angle in radians.
        translation: Translation vector (3,).

    Returns:
        Transformed node coordinates (N x 3).
    """
    if abs(angle) < 1e-12:
        return nodes + translation
    
    axis_normalized = axis / np.linalg.norm(axis)
    rotvec = axis_normalized * angle
    rotation = Rotation.from_rotvec(rotvec)
    rotated = rotation.apply(nodes)
    
    return rotated + translation