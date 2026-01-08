#ifndef PARTICLE_OPERATIONS_HPP
#define PARTICLE_OPERATIONS_HPP

#include <cmath>

#include "geometry.hpp"

namespace pygrain
{

constexpr double ROTATION_FACTOR = 1.0E-12;
constexpr double OVERLAP_TOLERANCE = 1.0E-6;

/**
 * @brief Translate a particle by a given vector.
 * 
 * @param geometry The packing geometry.
 * @param idx The particle index.
 * @param translation The translation vector [dx, dy, dz].
 */
void translate_particle(PackingGeometry& geometry,
                        std::size_t idx, 
                        const std::array<double, 3>& translation);

/**
 * @brief Set the absolute position of a particle.
 * 
 * @param geometry The packing geometry.
 * @param idx The particle index.
 * @param position The new position [x, y, z].
 */
void set_particle_position(PackingGeometry& geometry,
                           std::size_t idx,
                           const std::array<double, 3>& position);

/**
 * @brief Rotate a particle around a specified axis by a given angle.
 * 
 * @param geometry The packing geometry.
 * @param idx The particle index.
 * @param axis The rotation axis (assumed normalized) [ux, uy, uz].
 * @param angle The rotation angle in radians.
 */
void rotate_particle(PackingGeometry& geometry,
                     std::size_t idx, 
                     const std::array<double, 3>& axis, 
                     double angle);

/**
 * @brief Compute interaction between two particles and apply a displacement and rotation.
 * 
 * @param geometry The packing geometry.
 * @param idx1 The first particle index.
 * @param idx2 The second particle index.
 * @return Whether an overlap was detected.
 */
bool particle_interaction(PackingGeometry& geometry,
                            std::size_t idx1, 
                            std::size_t idx2);

} // namespace pygrain

#endif // PARTICLE_OPERATIONS_HPP