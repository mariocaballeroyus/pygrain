"""Surface roughness generation and application on meshes."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def assemble_autocorrelation_matrix(
    nodes: NDArray[np.float64], corr_length: float
) -> NDArray[np.float64]:
    """Assemble autocorrelation matrix for a set of nodes.
    
    Computes the autocorrelation matrix C where:
        C_{ij} = \exp{-d_{ij}^2 / (2 l^2)}
    
    where d_ij is the Euclidean distance between nodes i and j,
    and l is the correlation length.

    And computes its Cholesky decomposition C = L L^T.
    
    Args:
        nodes: Array of shape (N, 3) containing node coordinates.
        corr_length: Correlation length.
    
    Returns:
        Cholesky factor L of the autocorrelation matrix C.
    """
    d_ij = nodes[:, np.newaxis, :] - nodes[np.newaxis, :, :]
    d_ij_sq = np.sum(d_ij**2, axis=2)
    C = np.exp(-d_ij_sq / (2.0 * corr_length**2))

    try:
        L = np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        print("Warning: Autocorrelation matrix is not positive definite."
              " Applying regularization.")

        epsilon = 1e-10 * np.trace(C) / C.shape[0]
        C_reg = C + epsilon * np.eye(C.shape[0])
        
        try:
            L = np.linalg.cholesky(C_reg)
        except np.linalg.LinAlgError:
            raise ValueError("Failed to compute Cholesky decomposition")
    
    return L


def compute_roughness_profile(
    L: NDArray[np.float64], sq_roughness: float
) -> NDArray[np.float64]:
    """Generate a roughness profile from Cholesky factor.
    
    Args:
        L: Cholesky factor of autocorrelation matrix, shape (N, N).
        sq_roughness: Target squared roughness.
    
    Returns:
        Roughness profile of shape (N, 1).
    """
    # Gaussian white noise vector
    noise = np.random.normal(size=(L.shape[0], 1))

    roughness_profile = np.dot(L, noise)
    return roughness_profile * np.sqrt(sq_roughness)  # scaling


def compute_surface_normals(
    nodes: NDArray[np.float64], triangles: NDArray[np.int64]
) -> NDArray[np.float64]:
    """Compute surface normals at each node using area-weighted averaging.
    
    Args:
        nodes: Array of shape (N, 3) containing node coordinates.
        triangles: Array of shape (M, 3) containing triangle connectivity.
    
    Returns:
        Array of shape (N, 3) containing unit normal vectors at each node.
    """
    normals = np.zeros_like(nodes)
    
    # Get vertices for each triangle
    v0 = nodes[triangles[:, 0]]
    v1 = nodes[triangles[:, 1]]
    v2 = nodes[triangles[:, 2]]
    
    # Compute triangle normals (cross product of two edges)
    edge1 = v1 - v0
    edge2 = v2 - v0
    tri_normals = np.cross(edge1, edge2)
    
    # Accumulate normals to vertices (area-weighted)
    for i in range(3):
        np.add.at(normals, triangles[:, i], tri_normals)
    
    # Normalize to unit vectors
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    norms = np.where(norms > 1e-12, norms, 1.0)
    return normals / norms

def generate_and_apply_roughness(
    nodes: NDArray[np.float64],
    normals: NDArray[np.float64],
    L: NDArray[np.float64],
    sq_roughness: float
) -> NDArray[np.float64]:
    """Apply surface roughness using pre-computed Cholesky factor and normals.
    
    This function is optimized for cases where the Cholesky decomposition
    and surface normals have already been computed (e.g., for template meshes).
    
    Args:
        nodes: Array of shape (N, 3) containing node coordinates.
        normals: Array of shape (N, 3) containing unit normal vectors.
        cholesky_factor: Pre-computed Cholesky factor L of autocorrelation matrix.
        sq_roughness: Target squared roughness (variance).
    
    Returns:
        Modified nodes with roughness applied.
    """
    roughness_profile = compute_roughness_profile(L, sq_roughness)
    return nodes + roughness_profile * normals