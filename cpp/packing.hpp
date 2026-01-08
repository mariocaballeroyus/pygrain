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
    Packing() : geometry_(), lengths_({0.0, 0.0, 0.0}), geometry_idx_(0) {}

    /**
     * @brief Packing constructor with specified domain lengths.
     * 
     * @param lengths The domain box lengths in x, y, z directions.
     */
    Packing(std::array<double, 3> lengths)
        : geometry_(), lengths_(lengths), geometry_idx_(0) {}

    /**
     * @brief Add a spherical particle to the packing.
     * 
     * @param radius The sphere radius.
     * @param num Number of particles to add.
     * @param geometry_idx The geometry index to assign to these particles.
     */
    void add_sphere_particles(double radius, int num, std::size_t geometry_idx);

    /**
     * @brief Add a spheroidal particle to the packing.
     * @details Only prolate spheroids are supported (major axis >= minor axis).
     * 
     * @param aspect_ratio The aspect ratio (major/minor axis).
     * @param minor_axis The minor axis length.
     * @param num Number of particles to add.
     * @param geometry_idx The geometry index to assign to these particles.
     */
    void add_spheroid_particles(double aspect_ratio, double minor_axis, int num, std::size_t geometry_idx);

    /**
     * @brief Add a cylindrical particle to the packing.
     * @details Only prolate cylinders are supported (length >= diameter).
     * 
     * @param aspect_ratio The aspect ratio (length/diameter).
     * @param diameter The diameter.
     * @param num Number of particles to add.
     * @param geometry_idx The geometry index to assign to these particles.
     */
    void add_cylinder_particles(double aspect_ratio, double diameter, int num, std::size_t geometry_idx);

    /**
     * @brief Randomize positions and orientations of all particles.
     */
    void randomize_particles();

    /**
     * @brief Generate packing by resolving overlaps iteratively.
     * 
     * @param max_iterations Maximum number of iterations.
     */
    void generate(unsigned int max_iterations);

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
     * @brief Get the ID, position and orientation of a particle.
     * @details The orientation is expressed as axis-angle rotation from the 
     *          reference orientation (aligned with x-axis). This format is 
     *          compatible with Gmsh's rotation operations.
     * 
     * @param idx The particle index.
     * @return Array of 8 doubles: [id, x, y, z, axis_x, axis_y, axis_z, angle].
     */
    std::array<double, 8> particle_data(std::size_t idx) const;

    /**
     * @brief Get the positions and orientations of all particles.
     * @details Returns a flat array of size 8*N where N is the number of particles
     *          (or more if include_periodic is true).
     *          Each particle's data is [x, y, z, axis_x, axis_y, axis_z, angle, type].
     * 
     * @param include_periodic If true, includes periodic images that touch the primary box.
     * @return Vector of positions in row-major order.
     */
    std::vector<double> data_array(bool periodic) const;

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

    /**
     * @brief Counter for geometry group indices (increments with each add_*_particles call).
     */
    std::size_t geometry_idx_;
};

} // namespace pygrain

#endif // PACKING_HPP