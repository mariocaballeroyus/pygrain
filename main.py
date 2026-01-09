"""Example usage of pygrain with STL mesh export."""

import pygrain

def main():
    domain_size = 100.0
    packing = pygrain.packing.create_packing([domain_size, domain_size, domain_size])
    
    packing.add_sphere_particles(
        radius=11.5, corr_length=1.5, sq_roughness=0.1, num=30)

    packing.add_sphere_particles(
        radius=11.5, num=30)

    packing.add_spheroid_particles(
        aspect_ratio=2.0, minor_axis=10.0, corr_length=1.0, sq_roughness=0.1, num=35)

    packing.add_spheroid_particles(
        aspect_ratio=2.0, minor_axis=10.0, num=35)

    packing.add_cylinder_particles(
        aspect_ratio=2.0, diameter=15.0, num=30)

    print("\nGenerating packing...")
    packing.generate(max_iterations=50000)

    pygrain.io.to_csv(packing, "out/packing.csv")
    
    with pygrain.mesh.PackingMesh("packing") as mesh:
        mesh.generate(packing, mesh_size=0.7, periodic=True)
        pygrain.io.to_stl("out/packing_mesh.stl")

    print(f"\nPacking generated with:")
    print(f"  - {packing.num_particles} particles") 
    print(f"  - Porosity: {packing.porosity:.2%}")
    print(f"  - Surface area: {mesh.surface_area:.2f}")
    print(f"  - Specific surface: {mesh.surface_area / packing.volume:.2f}")


if __name__ == "__main__":
    main()
