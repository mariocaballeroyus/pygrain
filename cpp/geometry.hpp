#ifndef PARTICLE_GEOMETRY_HPP
#define PARTICLE_GEOMETRY_HPP

#include <vector>
#include <array>

namespace pygrain
{

/**
 * @brief Particle type enumeration.
 */
enum class ParticleType : int
{
    Sphere = 0,
    Spheroid = 1,
    Cylinder = 2
};

struct PackingGeometry
{
    /**
     * @brief Particle center of mass coordinates and bounding radii stored as contiguous arrays.
     * @details Format: [x, y, z, radius, idx]. 
     * All rotation operations are performed around a particle's center point. 
     * A bounding radius is defined as the radius of the smallest sphere enclosing the particle geometry.
     */
    std::array<std::vector<double>, 5> particle_data;

    /**
     * @brief Sphere center coordinates and radii from all particles stored as contiguous arrays.
     * @details Format: [x, y, z, radius].
     * Multiple spheres approximate the shape of each particle.
     */
    std::array<std::vector<double>, 4> sphere_data;

    /**
     * @brief Particle starting offsets in the sphere position and radii arrays.
     */
    std::vector<std::size_t> particle_offsets = {0};

    /** 
     * @brief Return the number of particles in the geometry.
     */
    std::size_t num_particles() const
    { return particle_data[0].size(); }

    /** 
     * @brief Return the total number of spheres in the geometry.
     */
    std::size_t num_spheres() const
    { return sphere_data[0].size(); }

    /** 
     * @brief Return the number of spheres comprising a given particle.
     * 
     * @param id The particle index.
     */
    std::size_t num_particle_spheres(std::size_t id) const
    { return particle_offsets[id + 1] - particle_offsets[id]; }
};

} // namespace pygrain

#endif // PARTICLE_GEOMETRY_HPP