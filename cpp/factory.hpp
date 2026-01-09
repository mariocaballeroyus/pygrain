#ifndef FACTORY_HPP
#define FACTORY_HPP

#include "geometry.hpp"

namespace pygrain3d
{

/**
 * @brief Default sphere placement precision for particle approximations.
 * @details Defines the density of spheres used to approximate complex particle shapes.
 *          Higher values yield better approximations at a computational cost.
 */
constexpr double SPHERE_PRECISION = 2.0;

/**
 * @brief Begin defining a new particle in the geometry.
 * 
 * @param geometry The packing geometry.
 * @param px Particle center x-coordinate.
 * @param py Particle center y-coordinate.
 * @param pz Particle center z-coordinate.
 * @param pr Particle bounding radius.
 * @param id Geometry group ID (particles with same geometry share this).
 */
void begin_particle(PackingGeometry& geometry, 
                    double px, double py, double pz, 
                    double pr, 
                    std::size_t id);

/** 
 * @brief Finalize the current particle definition.
 * @details Records the current number of spheres as the end offset for the particle.
 * 
 * @param geometry The packing geometry.
 */
void end_particle(PackingGeometry& geometry);

/**
 * @brief Add a sphere to the current particle definition.
 * @details The sphere is centered at the origin.
 * 
 * @param geometry The packing geometry.
 * @param radius Sphere radius.
 * @param id Geometry group ID (particles with same geometry share this).
 */
void generate_sphere_particle(PackingGeometry& geometry, 
                              double radius,
                              std::size_t id);

/**
 * @brief Generate a spheroidal particle approximated by spheres.
 * @details Only prolate spheroids are supported (major axis >= minor axis).
 *          The spheroid is centered at the origin and aligned along the x-axis.
 * 
 * @param geometry The packing geometry.
 * @param aspect_ratio The aspect ratio (major/minor axis).
 * @param minor_diameter The minor axis length.
 * @param id Geometry group ID (particles with same geometry share this).
 */
void generate_spheroid_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double minor_diameter,
                                std::size_t id);

/**
 * @brief Generate a cylindrical particle approximated by spheres.
 * @details Only prolate cylinders are supported (length >= diameter).
 *          The cylinder is centered at the origin and aligned along the x-axis.
 * 
 * @param geometry The packing geometry.
 * @param aspect_ratio The aspect ratio (length/diameter).
 * @param diameter The diameter.
 * @param id Geometry group ID (particles with same geometry share this).
 */
void generate_cylinder_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double diameter,
                                std::size_t id);

} // namespace pygrain

#endif // FACTORY_HPP