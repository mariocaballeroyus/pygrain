#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <array>

namespace pygrain3d
{

struct PackingGeometry
{
    /**
     * @brief Particle ID, center of mass coordinates and bounding radii stored as contiguous arrays.
     * @details Format: [id, x, y, z, r]. 
     *          All rotation operations are performed around a particle's center point. 
     *          A bounding radius is defined as the radius of the smallest sphere enclosing the 
     *          particle geometry.
     */
    std::array<std::vector<double>, 5> particles;

    /**
     * @brief Sphere center coordinates and radii from all particles stored as contiguous arrays.
     * @details Format: [x, y, z, radius].
     *          Multiple spheres approximate the shape of each particle.
     */
    std::array<std::vector<double>, 4> spheres;

    /**
     * @brief Particle starting offsets in the sphere position and radii arrays.
     */
    std::vector<std::size_t> particle_offsets = {0};

    /** 
     * @brief Return the number of particles in the geometry.
     */
    std::size_t num_particles() const
    { return particles[0].size(); }

    /** 
     * @brief Return the total number of spheres in the geometry.
     */
    std::size_t num_spheres() const
    { return spheres[0].size(); }

    /** 
     * @brief Return the number of spheres comprising a given particle.
     * 
     * @param id The particle index.
     */
    std::size_t num_particle_spheres(std::size_t id) const
    { return particle_offsets[id + 1] - particle_offsets[id]; }

    /**
     * @brief Get the ID, position and orientation of a particle.
     * @details A particle's data is [id, x, y, z, axis_x, axis_y, axis_z, angle].
     * 
     * @param idx The particle index.
     * @return The array of particle data.
     */
    std::array<double, 8> particle_data(std::size_t idx) const;

    /**
     * @brief Get the positions and orientations of all particles.
     * @details Each particle's data is [id, x, y, z, axis_x, axis_y, axis_z, angle].
     * 
     * @param lengths The domain box lengths for periodic image calculation.
     * @param periodic Whether to include periodic images potentially overlapping the box.
     * @return The array of all particle data.
     */
    std::vector<double> particle_array(const std::array<double, 3>& lengths, bool periodic) const;

    /**
     * @brief Get the positions and radii of all spheres.
     * @details Each sphere's data is [x, y, z, radius].
     * 
     * @return The array of all sphere data.
     */
    std::vector<double> sphere_array() const;

    /**
     * @brief Clear all particles and spheres from the geometry.
     * @details Removes all particle and sphere data, resetting the geometry to empty state.
     */
    void clear();
};

} // namespace pygrain3d

#endif // GEOMETRY_HPP