#include "operations.hpp"

namespace pygrain
{

void translate_particle(PackingGeometry& geometry,
                        std::size_t idx, 
                        const std::array<double, 3>& translation)
{
    auto& [px, py, pz, pr, pid] = geometry.particle_data;
    auto& [sx, sy, sz, sr] = geometry.sphere_data;

    std::size_t start = geometry.particle_offsets[idx];
    std::size_t end = geometry.particle_offsets[idx + 1];

    // Update particle center position
    px[idx] += translation[0];
    py[idx] += translation[1];
    pz[idx] += translation[2];

    // Update sphere positions
    for (std::size_t i = start; i < end; ++i) 
    {
        sx[i] += translation[0];
        sy[i] += translation[1];
        sz[i] += translation[2];
    }
}

void set_particle_position(PackingGeometry& geometry,
                           std::size_t idx,
                           const std::array<double, 3>& position)
{
    auto& [px, py, pz, pr, pid] = geometry.particle_data;
    auto& [sx, sy, sz, sr] = geometry.sphere_data;

    std::size_t start = geometry.particle_offsets[idx];
    std::size_t end = geometry.particle_offsets[idx + 1];

    // Compute translation vector
    const double dx = position[0] - px[idx];
    const double dy = position[1] - py[idx];
    const double dz = position[2] - pz[idx];

    // Update particle center position
    px[idx] = position[0];
    py[idx] = position[1];
    pz[idx] = position[2];

    // Update sphere positions
    for (std::size_t i = start; i < end; ++i) 
    {
        sx[i] += dx;
        sy[i] += dy;
        sz[i] += dz;
    }
}

void rotate_particle(PackingGeometry& geometry,
                     std::size_t idx, 
                     const std::array<double, 3>& axis, 
                     double angle)
{
    auto& [px, py, pz, pr, pid] = geometry.particle_data;
    auto& [sx, sy, sz, sr] = geometry.sphere_data;
    const auto& [ux, uy, uz] = axis;  // assumed normalized

    // Pre-calculate trigonometry
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double one_c = 1.0 - c;

    // Pre-calculate rotation matrix
    const double r00 = c + ux * ux * one_c;         // Row 1
    const double r01 = ux * uy * one_c - uz * s;
    const double r02 = ux * uz * one_c + uy * s;
    const double r10 = uy * ux * one_c + uz * s;    // Row 2
    const double r11 = c + uy * uy * one_c;
    const double r12 = uy * uz * one_c - ux * s;
    const double r20 = uz * ux * one_c - uy * s;    // Row 3
    const double r21 = uz * uy * one_c + ux * s;
    const double r22 = c + uz * uz * one_c;

    const std::size_t start = geometry.particle_offsets[idx];
    const std::size_t end = geometry.particle_offsets[idx + 1];

    const double cx = px[idx];
    const double cy = py[idx];
    const double cz = pz[idx];

    for (std::size_t i = start; i < end; ++i) 
    {
        // Translate to origin
        const double x = sx[i] - cx;
        const double y = sy[i] - cy;
        const double z = sz[i] - cz;

        // Apply rotation matrix
        const double rx = r00 * x + r01 * y + r02 * z;
        const double ry = r10 * x + r11 * y + r12 * z;
        const double rz = r20 * x + r21 * y + r22 * z;

        // Translate back
        sx[i] = rx + cx;
        sy[i] = ry + cy;
        sz[i] = rz + cz;
    }
}

bool particle_interaction(PackingGeometry& geometry,
                            std::size_t idx1, 
                            std::size_t idx2)
{
    auto& [px, py, pz, pr, pid] = geometry.particle_data;
    auto& [sx, sy, sz, sr] = geometry.sphere_data;

    std::array<double, 3> displacement = {0.0, 0.0, 0.0};
    std::array<double, 3> rotation = {0.0, 0.0, 0.0};
    double max_overlap_mag = 0.0;
    std::array<double, 3> max_overlap_dir = {0.0, 0.0, 0.0};

    const std::size_t start1 = geometry.particle_offsets[idx1];
    const std::size_t end1 = geometry.particle_offsets[idx1 + 1];
    const std::size_t start2 = geometry.particle_offsets[idx2];
    const std::size_t end2 = geometry.particle_offsets[idx2 + 1];

    for (std::size_t m = start1; m < end1; ++m)  // spheres of particle 1
    {
        for (std::size_t n = start2; n < end2; ++n)  // spheres of particle 2
        {
            // Compute overlap vector between spheres m and n
            const double dx = sx[n] - sx[m];
            const double dy = sy[n] - sy[m];
            const double dz = sz[n] - sz[m];
            const double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            const double radius_sum = sr[m] + sr[n];
            const double tol = OVERLAP_TOLERANCE * radius_sum;
            
            // Only count as overlap if exceeds tolerance
            if (dist < radius_sum - tol)
            {
                const double overlap_mag = radius_sum - dist;
                
                // Keep track of the maximum overlap
                if (overlap_mag > max_overlap_mag)
                {
                    max_overlap_mag = overlap_mag;
                    // Direction points from n to m (push m away from n)
                    max_overlap_dir[0] = -dx / dist;
                    max_overlap_dir[1] = -dy / dist;
                    max_overlap_dir[2] = -dz / dist;
                }
            }
        }
    }

    if (max_overlap_mag == 0.0)
        return false;  // no overlap detected

    // Displacement based on maximum overlap
    displacement[0] = max_overlap_dir[0] * max_overlap_mag;
    displacement[1] = max_overlap_dir[1] * max_overlap_mag;
    displacement[2] = max_overlap_dir[2] * max_overlap_mag;

    const double displacement_magnitude = max_overlap_mag;

    translate_particle(geometry, idx1, displacement);

    // Rotation: cross product (center_i - center_j) x displacement
    const double Ci_x = px[idx1];
    const double Ci_y = py[idx1];
    const double Ci_z = pz[idx1];
    const double Cj_x = px[idx2];
    const double Cj_y = py[idx2];
    const double Cj_z = pz[idx2];

    // Compute torque vector
    rotation[0] = (Ci_y - Cj_y) * displacement[2] - (Ci_z - Cj_z) * displacement[1];
    rotation[1] = (Ci_z - Cj_z) * displacement[0] - (Ci_x - Cj_x) * displacement[2];
    rotation[2] = (Ci_x - Cj_x) * displacement[1] - (Ci_y - Cj_y) * displacement[0];

    // Compute rotation angle (non-scaled)
    const double rot_angle = std::sqrt(rotation[0]*rotation[0] +
                                       rotation[1]*rotation[1] +
                                       rotation[2]*rotation[2]);

    if (rot_angle > 1e-8)
    {
        std::array<double, 3> rot_axis = {rotation[0]/rot_angle,
                                          rotation[1]/rot_angle,
                                          rotation[2]/rot_angle};

        rotate_particle(geometry, idx1, rot_axis, rot_angle * ROTATION_FACTOR);
    }

    return true;
}

} // namespace pygrain