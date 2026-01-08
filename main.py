#!/usr/bin/env python3
"""Example usage of pygrain with STL mesh export."""

import pygrain

def main():
    # Create a packing domain
    domain_size = 100.0
    packing = pygrain.packing.create_packing([domain_size, domain_size, domain_size])
    
    # Add particles in batches - each batch shares the same geometry
    packing.add_sphere_particles(radius=11.5, num=60)
    packing.add_spheroid_particles(aspect_ratio=2.0, minor_axis=10.0, num=70)
    packing.add_cylinder_particles(aspect_ratio=2.0, diameter=15.0, num=30)

    print("\nGenerating packing...")
    packing.generate(max_iterations=50000)

    print(f"Packing generated with:")
    print(f"\t{packing.num_particles} particles") 
    print(f"\tPorosity: {packing.porosity:.2%}")

    pygrain.io.to_csv(packing, "out/packing.csv")
    
    with pygrain.Mesh("packing") as mesh:
        mesh.generate_from_packing(packing, mesh_size=2.0, periodic=True)
        pygrain.io.to_stl(mesh, "out/packing_mesh.stl")


if __name__ == "__main__":
    main()
