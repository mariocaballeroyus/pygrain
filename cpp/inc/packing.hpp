#ifndef PACKING_HPP
#define PACKING_HPP

#include "packing_geometry.hpp"
#include "particle_factory.hpp"

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
     */
    void add_sphere_particle(double radius);

    /**
     * @brief Add a spheroidal particle to the packing.
     * @details Only prolate spheroids are supported (major axis >= minor axis).
     * 
     * @param aspect_ratio The aspect ratio (major/minor axis).
     * @param minor_axis The minor axis length.
     */
    void add_spheroid_particle(double aspect_ratio, double minor_axis);

    /**
     * @brief Add a cylindrical particle to the packing.
     * @details Only prolate cylinders are supported (length >= diameter).
     * 
     * @param aspect_ratio The aspect ratio (length/diameter).
     * @param diameter The diameter.
     */
    void add_cylinder_particle(double aspect_ratio, double diameter);

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
     * @brief Export spheres to CSV file in x, y, z, R format.
     * 
     * @param filename The output filename.
     */
    void export_spheres_csv(const std::string& filename) const;

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