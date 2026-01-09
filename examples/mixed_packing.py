# %% [markdown]
# ## Mixed Particle Packing

# %%
import pygrain3d as pg
print(pg.__version__)

# %%
length = 100.0
packing = pg.packing.create_packing([3*length, length, length])

# %%
packing.add_sphere_particles(
    radius=11.5, corr_length=1.0, sq_roughness=0.1, num=45)

packing.add_sphere_particles(
    radius=11.5, corr_length=2.0, sq_roughness=0.3, num=45)

packing.add_sphere_particles(
    radius=11.5, num=90)

packing.add_spheroid_particles(
    aspect_ratio=2.0, minor_axis=10.0, corr_length=1.0, sq_roughness=0.1, num=50)

packing.add_spheroid_particles(
    aspect_ratio=2.0, minor_axis=10.0, corr_length=2.0, sq_roughness=0.2, num=50)

packing.add_spheroid_particles(
    aspect_ratio=2.0, minor_axis=10.0, num=105)

packing.add_cylinder_particles(
    aspect_ratio=2.0, diameter=15.0, num=90)

# %%
packing.generate(max_iterations=50000)
pg.io.to_csv(packing, "mixed_packing.csv")

# %%
with pg.mesh.PackingMesh("packing") as mesh:
    mesh.generate(packing, mesh_size=0.7, periodic=True)
    pg.io.to_stl("mixed_packing_mesh.stl")