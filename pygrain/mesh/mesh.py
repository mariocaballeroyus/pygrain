"""Packing mesh generation using Gmsh."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

import gmsh

from ..packing.packing import Packing
from ..packing.particle import Sphere, Spheroid
from .geometry import (
    _create_occ_geometry, 
    _mesh_occ_geometry, 
    _transform_nodes
)
from .roughness import (
    assemble_autocorrelation_matrix,
    compute_surface_normals,
    generate_and_apply_roughness
)


class PackingMesh:
    """Packing surface mesh generator."""

    def __init__(self, name: str = "model") -> None:
        """Initialize a `PackingMesh` instance.
        
        Args:
            name: Name of the Gmsh model (default: "model").
        """
        self._name = name
        self._nodes: NDArray[np.float64] | None = None
        self._triangles: NDArray[np.int64] | None = None

        self._num_particles: int = 0  # track particles without periodic images
        self._num_triangles: int = 0  # non-image triangles to compute area
        self._surface_area: float = 0.0
    
    def __enter__(self) -> "PackingMesh":
        """Enter context manager, initialize Gmsh."""
        gmsh.initialize()
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.model.add(self._name)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Exit context manager, finalize Gmsh."""
        gmsh.finalize()
        return False
    
    def _check_gmsh_initialized(self) -> None:
        """Check if gmsh is initialized, raise error if not."""
        if not gmsh.isInitialized():
            raise RuntimeError(
                "Gmsh is not initialized. "
                "Use `PackingMesh` as a context manager or call gmsh.initialize()."
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
        
        self._nodes, self._triangles, self._num_triangles = self._generate_combined_mesh(
            packing, mesh_size, periodic=periodic
        )
        
        # Compute surface area
        self._surface_area = self._compute_surface_area()
        
        self._populate_gmsh_model()
        
        print(f"Generated mesh with {len(self._nodes)} nodes, {len(self._triangles)} triangles")
        print(f"Total surface area (primary particles): {self._surface_area:.2f}")
        
    def _generate_combined_mesh(
        self, packing: Packing, mesh_size: float, periodic: bool = False
    ) -> tuple[NDArray[np.float64], NDArray[np.int64], int, int]:
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
            Tuple of (all_nodes, all_triangles, num_primary_particles, num_primary_triangles).
        """
        data_array = packing.data_array(periodic=periodic)
        
        if len(data_array) == 0:
            raise ValueError("No particles in packing")
        
        # Create one template mesh per unique geometry. Stores:
        # - Nodes
        # - Triangles
        # - Cholesky factor (or None if no roughness)
        # - Normals (or None if no roughness)
        templates: dict[int, tuple[NDArray[np.float64], 
                                   NDArray[np.int64], 
                                   NDArray[np.float64] | None, 
                                   NDArray[np.float64] | None]] = {}
        
        for idx, particle in enumerate(packing.particles):
            gmsh.model.add(f"template_{idx}")
            
            vol_tag = _create_occ_geometry(particle)
            nodes, triangles = _mesh_occ_geometry(mesh_size)
            
            # Compute Cholesky and normals if particle has roughness
            cholesky_factor = None
            template_normals = None
            if isinstance(particle, (Sphere, Spheroid)):
                if particle.corr_length > 0 and particle.sq_roughness > 0:
                    cholesky_factor = assemble_autocorrelation_matrix(nodes, particle.corr_length)
                    template_normals = compute_surface_normals(nodes, triangles)
            
            templates[idx] = (nodes, triangles, cholesky_factor, template_normals)
            gmsh.model.remove()
        
        # Duplicate and transform for each particle instance
        all_nodes_list: list[NDArray[np.float64]] = []
        all_elements_list: list[NDArray[np.int64]] = []
        node_offset = 0
        triangle_offset = 0
        
        # Track which particles we've seen (for periodic image handling)
        # Periodic images have same position in primary box, so they share roughness
        roughness_cache: dict[tuple[int, int], NDArray[np.float64]] = {}
        
        # Track number of triangles for primary particles
        num_primary_triangles = 0
        particle_count = 0
        
        for data in data_array:
            geometry_idx = int(data[0])
            translation = data[1:4]
            axis = data[4:7]
            angle = data[7]
            
            template_nodes, template_elements, cholesky_factor, template_normals = templates[geometry_idx]
            
            # Apply roughness if needed
            if cholesky_factor is not None and template_normals is not None:
                # Periodic images will have different positions but identical mesh
                # We need to map back to primary box position
                primary_pos = translation % np.array([packing.lengths[0], packing.lengths[1], packing.lengths[2]])
                pos_hash = hash(tuple(np.round(primary_pos, decimals=6)))
                cache_key = (geometry_idx, pos_hash)
                
                if cache_key in roughness_cache:
                    # Reuse roughness for periodic images
                    nodes_to_transform = roughness_cache[cache_key]
                else:
                    # Generate unique roughness for this particle
                    particle = packing.particles[geometry_idx]
                    
                    nodes_to_transform = generate_and_apply_roughness(
                        template_nodes, 
                        template_normals, 
                        cholesky_factor, 
                        particle.sq_roughness
                    )
                    
                    # Cache for periodic images
                    roughness_cache[cache_key] = nodes_to_transform
            else:
                nodes_to_transform = template_nodes
            
            transformed_nodes = _transform_nodes(nodes_to_transform, axis, angle, translation)
            offset_elements = template_elements + node_offset
            
            all_nodes_list.append(transformed_nodes)
            all_elements_list.append(offset_elements)
            
            # Track primary particles (count up to num_particles)
            if particle_count < packing.num_particles:
                num_primary_triangles += len(template_elements)
                particle_count += 1
            
            node_offset += len(transformed_nodes)
            triangle_offset += len(template_elements)
        
        all_nodes = np.vstack(all_nodes_list)
        all_triangles = np.vstack(all_elements_list)
                
        return all_nodes, all_triangles, num_primary_triangles
    
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

    def _compute_surface_area(self) -> float:
        """Compute total surface area of primary particles (excluding periodic images).
        
        Returns:
            Total surface area.
        """
        if self._nodes is None or self._triangles is None or self._num_triangles == 0:
            return 0.0
        
        # Only compute surface area for primary particles (first num_triangles)
        primary_triangles = self._triangles[:self._num_triangles]
        
        # Get triangle vertices
        v0 = self._nodes[primary_triangles[:, 0]]
        v1 = self._nodes[primary_triangles[:, 1]]
        v2 = self._nodes[primary_triangles[:, 2]]
        
        # Compute triangle areas using cross product
        edge1 = v1 - v0
        edge2 = v2 - v0
        cross = np.cross(edge1, edge2)
        areas = 0.5 * np.linalg.norm(cross, axis=1)
        
        # Sum all triangle areas
        total_area = np.sum(areas)
        
        return total_area
    
    @property
    def surface_area(self) -> float:
        """Total surface area of primary particles (excluding periodic images)."""
        return self._surface_area