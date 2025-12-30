#ifndef PARTICLE_FACTORY_HPP
#define PARTICLE_FACTORY_HPP

#include "packing_geometry.hpp"

namespace pygrain
{

/**
 * @brief Begin defining a new particle in the geometry.
 * 
 * @param geometry The packing geometry.
 * @param px Particle center x-coordinate.
 * @param py Particle center y-coordinate.
 * @param pz Particle center z-coordinate.
 * @param bounding_radius Particle bounding radius.
 */
void begin_particle(PackingGeometry& geometry, 
                    double px, double py, double pz,
                    double bounding_radius);

/** 
 * @brief Finalize the current particle definition.
 * @details Records the current number of spheres as the end offset for the particle.
 * 
 * @param geometry The packing geometry.
 */
void end_particle(PackingGeometry& geometry);

/**
 * @brief Add a sphere to the current particle definition.
 * @details The sphere particle position is set at the origin.
 * 
 * @param geometry The packing geometry.
 * @param radius Sphere radius.
 */
void generate_sphere_particle(PackingGeometry& geometry, 
                              double radius);

/**
 * @brief Generate a spheroidal particle approximated by spheres.
 * @details Only prolate spheroids are supported (major axis >= minor axis).
 *          The spheroid is centered at the origin and aligned along the x-axis.
 * 
 * @param geometry The packing geometry.
 * @param aspect_ratio The aspect ratio (major/minor axis).
 * @param minor_diameter The minor axis length.
 * @param sphere_precision Controls the density of spheres along the major axis.
 */
void generate_spheroid_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double minor_diameter, 
                                double sphere_precision = 2.0);

/**
 * @brief Generate a cylindrical particle approximated by spheres.
 * @details Only prolate cylinders are supported (length >= diameter).
 *          The cylinder is centered at the origin and aligned along the x-axis.
 * 
 * @param geometry The packing geometry.
 * @param aspect_ratio The aspect ratio (length/diameter).
 * @param diameter The diameter.
 * @param sphere_precision Controls the density of spheres along the cylinder axis.
 */
void generate_cylinder_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double diameter, 
                                double sphere_precision = 2.0);

} // namespace pygrain

#endif // PARTICLE_FACTORY_HPP