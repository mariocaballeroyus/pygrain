#ifndef PACKING_HPP
#define PACKING_HPP

#include "geometry.hpp"
#include "factory.hpp"

#include <memory>
#include <array>
#include <random>
#include <iostream>
#include <fstream>

namespace pygrain
{

class Packing
{
public:
    /**
     * @brief Default constructor.
     * @details Initializes an empty packing with zero domain lengths.
     */
    Packing() : geometry_(), lengths_({0.0, 0.0, 0.0}) {}

    /**
     * @brief Packing constructor with specified domain lengths.
     * 
     * @param lengths The domain box lengths in x, y, z directions.
     */
    Packing(std::array<double, 3> lengths)
        : geometry_(), lengths_(lengths) {}

    /**
     * @brief Add a spherical particle to the packing.
     * 
     * @param radius The sphere radius.
     * @param num Number of particles to add.
     * @param id Geometry group ID (particles with same geometry share this).
     */
    void add_sphere_particles(double radius, 
                              int num, 
                              std::size_t id);

    /**
     * @brief Add a spheroidal particle to the packing.
     * @details Only prolate spheroids are supported (major axis >= minor axis).
     * 
     * @param aspect_ratio The aspect ratio (major/minor axis).
     * @param minor_axis The minor axis length.
     * @param num Number of particles to add.
     * @param id Geometry group ID (particles with same geometry share this).
     */
    void add_spheroid_particles(double aspect_ratio, 
                                double minor_axis, 
                                int num, 
                                std::size_t id);

    /**
     * @brief Add a cylindrical particle to the packing.
     * @details Only prolate cylinders are supported (length >= diameter).
     * 
     * @param aspect_ratio The aspect ratio (length/diameter).
     * @param diameter The diameter.
     * @param num Number of particles to add.
     * @param id Geometry group ID (particles with same geometry share this).
     */
    void add_cylinder_particles(double aspect_ratio, 
                                double diameter, 
                                int num, 
                                std::size_t id);

    /**
     * @brief Randomize positions and orientations of all particles in the packing.
     */
    void randomize_particles();

    /**
     * @brief Generate packing by resolving overlaps iteratively.
     * 
     * @param max_iterations Maximum number of iterations.
     */
    void generate(unsigned int max_iterations,
                  unsigned int log_interval);

    /**
     * @brief Get the ID, position and orientation of a particle.
     * @details A particle's data is [x, y, z, axis_x, axis_y, axis_z, angle, type].
     * 
     * @param idx The particle index.
     * @return The array of particle data.
     */
    std::array<double, 8> particle_data(std::size_t idx) const;

    /**
     * @brief Get the positions and orientations of all particles.
     * @details Each particle's data is [id, x, y, z, axis_x, axis_y, axis_z, angle].
     * 
     * @param periodic Whether to include periodic images potentially overlapping the box.
     * @return The array of all particle data.
     */
    std::vector<double> data_array(bool periodic) const;

    /**
     * @brief Get the (const) geometry of the packing.
     */
    const PackingGeometry& get_geometry() const 
    { return geometry_; }

    /**
     * @brief Get the (non-const) geometry of the packing.
     */
    PackingGeometry& get_geometry() 
    { return geometry_; }

    /**
     * @brief Get the domain lengths.
     */
    const std::array<double, 3>& get_lengths() const 
    { return lengths_; }

    /**
     * @brief Get the number of particles in the packing.
     */
    std::size_t num_particles() const
    { return geometry_.num_particles(); }

private:
    /**
     * @brief The packing geometry containing particles and spheres.
     */
    PackingGeometry geometry_;

    /**
     * @brief The domain box lengths in x, y, z directions.
     */
    std::array<double, 3> lengths_;
};

} // namespace pygrain

#endif // PACKING_HPP