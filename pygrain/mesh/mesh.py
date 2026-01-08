from __future__ import annotations

from contextlib import contextmanager
from typing import Generator

import numpy as np
from numpy.typing import NDArray

import gmsh

from ..packing.packing import Packing
from .geometry import (
    _create_occ_geometry, 
    _mesh_occ_geometry, 
    _transform_nodes
)


class PackingMesh:

    def __init__(self, name: str = "model") -> None:
        self._name = name
        self._nodes: NDArray[np.float64] | None = None
        self._triangles: NDArray[np.int64] | None = None
    
    def __enter__(self) -> "PackingMesh":
        gmsh.initialize()
        gmsh.model.add(self._name)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        gmsh.finalize()
        return False
    
    def _check_gmsh_initialized(self) -> None:
        """Check if gmsh is initialized, raise error if not."""
        if not gmsh.isInitialized():
            raise RuntimeError(
                "Gmsh is not initialized. "
                "Use `PackingMesh` as a context manager or call gmsh.initialize()"
            )

    def generate(
        self, packing: Packing, mesh_size: float, periodic: bool = True
    ) -> None:
        """Generate a surface mesh of the particle packing.
        
        Args:
            packing: The packing object containing particles.
            mesh_size: Target mesh element size.
            periodic: If True, includes periodic images that touch the box.
        """
        self._check_gmsh_initialized()
        
        self._nodes, self._triangles = self._generate_combined_mesh(
            packing, mesh_size, periodic=periodic
        )
        
        self._populate_gmsh_model()
        
        print(f"Generated mesh with {len(self._nodes)} nodes, {len(self._triangles)} triangles")
        
    def _generate_combined_mesh(
        packing: "Packing",
        mesh_size: float,
        periodic: bool = False
    ) -> tuple[NDArray[np.float64], NDArray[np.int64]]:
        """Generate combined mesh by creating templates and duplicating them.
        
        Workflow:
        1. Create OCC geometry for each particle type
        2. Mesh each template
        3. Store nodes and elements
        4. Duplicate and transform for each particle instance
        
        Args:
            packing: The packing object containing particles and positions.
            mesh_size: Target mesh element size.
            periodic: If True, includes periodic images that touch the box.
        
        Returns:
            Tuple of (all_nodes, all_triangles) for the combined mesh.
        """
        data_array = packing.data_array(periodic=periodic)
        
        if len(data_array) == 0:
            raise ValueError("No particles in packing")
        
        # Create one template mesh per unique geometry
        templates: dict[int, tuple[NDArray[np.float64], NDArray[np.int64]]] = {}
        
        for idx, particle in enumerate(packing.particles):
            gmsh.model.add(f"template_{idx}")
            
            vol_tag = _create_occ_geometry(particle)
            nodes, triangles = _mesh_occ_geometry(vol_tag, mesh_size)
            
            templates[idx] = (nodes, triangles)
            gmsh.model.remove()
        
        # Duplicate and transform for each particle instance
        all_nodes_list: list[NDArray[np.float64]] = []
        all_elements_list: list[NDArray[np.int64]] = []
        node_offset = 0
        
        for data in data_array:
            geometry_idx = int(data[0])
            translation = data[1:4]
            axis = data[4:7]
            angle = data[7]
            
            template_nodes, template_elements = templates[geometry_idx]
            
            transformed_nodes = _transform_nodes(template_nodes, axis, angle, translation)
            offset_elements = template_elements + node_offset
            
            all_nodes_list.append(transformed_nodes)
            all_elements_list.append(offset_elements)
            node_offset += len(transformed_nodes)
        
        all_nodes = np.vstack(all_nodes_list)
        all_triangles = np.vstack(all_elements_list)
        
        return all_nodes, all_triangles

    def _populate_gmsh_model(self) -> None:
        """Populate the Gmsh model with the generated nodes and triangles."""
        if self._nodes is None or self._triangles is None:
            return

        try:
            gmsh.model.setCurrent(self._name)
        except Exception:
            gmsh.model.add(self._name)
            
        gmsh.model.remove()
        gmsh.model.add(self._name)

        s_tag = gmsh.model.addDiscreteEntity(2)

        num_nodes = len(self._nodes)
        node_tags = list(range(1, num_nodes + 1))
        flat_coords = self._nodes.flatten()
        
        gmsh.model.mesh.addNodes(2, s_tag, node_tags, flat_coords)
        
        if len(self._triangles) > 0:
            elem_type = 2
            num_tris = len(self._triangles)
            elem_tags = list(range(1, num_tris + 1))
            flat_elem_node_tags = (self._triangles + 1).flatten()
            
            gmsh.model.mesh.addElements(2, s_tag, [elem_type], [elem_tags], [flat_elem_node_tags])