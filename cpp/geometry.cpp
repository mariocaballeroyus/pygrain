#include "geometry.hpp"

#include <algorithm>
#include <cmath>

namespace pygrain3d
{

std::array<double, 8> PackingGeometry::particle_data(std::size_t idx) const
{
    const auto& [px, py, pz, pr, pid] = particles;
    const auto& [sx, sy, sz, sr] = spheres;
    const auto& offsets = particle_offsets;

    int id = static_cast<int>(pid[idx]);
    double x = px[idx];
    double y = py[idx];
    double z = pz[idx];

    const std::size_t start = offsets[idx];
    const std::size_t end = offsets[idx + 1];
    
    double axis_x = 0.0, axis_y = 0.0, axis_z = 1.0;
    double angle = 0.0;

    if ((end - start) >= 2)  // two spheres needed to define orientation
    {
        // Define direction as the vector from sphere 0 to sphere 1
        double dir_x = sx[start + 1] - sx[start];
        double dir_y = sy[start + 1] - sy[start];
        double dir_z = sz[start + 1] - sz[start];

        double norm = std::sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
        if (norm > 1e-12)
        {
            dir_x /= norm;
            dir_y /= norm;
            dir_z /= norm;

            // Compute rotation from reference (1,0,0) to current direction
            double dot = std::max(-1.0, std::min(1.0, dir_x));
            angle = std::acos(dot);

            // Cross product: (1,0,0) x (dir_x, dir_y, dir_z)
            axis_x = 0.0;
            axis_y = -dir_z;
            axis_z = dir_y;

            double axis_norm = std::sqrt(axis_y*axis_y + axis_z*axis_z);
            if (axis_norm > 1e-12)
            {
                axis_y /= axis_norm;
                axis_z /= axis_norm;
            }
            else if (dot < 0)  // anti-parallel: direction is (-1,0,0)
            {
                axis_y = 1.0;
                axis_z = 0.0;
                angle = M_PI;
            }
        }
    }
    else  // no orientation, default to zero rotation
    {
        axis_x = 0.0;
        axis_y = 0.0;
        axis_z = 1.0;
        angle = 0.0;
    }

    return {static_cast<double>(id), x, y, z, axis_x, axis_y, axis_z, angle};
}

std::vector<double> PackingGeometry::particle_array(const std::array<double, 3>& lengths, bool periodic) const
{
    const std::size_t n = num_particles();
    std::vector<double> positions;
    
    // Optimization: Reserve space for at least the base particles
    positions.reserve(n * 8);

    // Add all base particles regardless of periodic value
    for (std::size_t i = 0; i < n; ++i)
    {
        auto pos = particle_data(i); 
        for (std::size_t j = 0; j < 8; ++j)
        {
            positions.push_back(pos[j]);
        }
    }

    // Handle Periodic Images
    if (periodic)
    {
        const auto& [px, py, pz, pr, pid] = particles;

        for (std::size_t i = 0; i < n; ++i)
        {
            const double radius = pr[i];
            
            // Determine possible shifts for each dimension
            // A shift is needed if the center is within 'radius' of a boundary
            std::vector<double> dx_shifts = {0.0};
            if (px[i] < radius) dx_shifts.push_back(lengths[0]);
            else if (px[i] > lengths[0] - radius) dx_shifts.push_back(-lengths[0]);

            std::vector<double> dy_shifts = {0.0};
            if (py[i] < radius) dy_shifts.push_back(lengths[1]);
            else if (py[i] > lengths[1] - radius) dy_shifts.push_back(-lengths[1]);

            std::vector<double> dz_shifts = {0.0};
            if (pz[i] < radius) dz_shifts.push_back(lengths[2]);
            else if (pz[i] > lengths[2] - radius) dz_shifts.push_back(-lengths[2]);

            // If we only have {0,0,0}, no images are needed
            if (dx_shifts.size() == 1 && dy_shifts.size() == 1 && dz_shifts.size() == 1)
                continue;

            // Cache orientation/ID to avoid repeated function calls
            auto base_data = particle_data(i);

            // Loop through combinations to generate 1, 3, or 7 images
            for (double dx : dx_shifts)
            {
                for (double dy : dy_shifts)
                {
                    for (double dz : dz_shifts)
                    {
                        if (dx == 0.0 && dy == 0.0 && dz == 0.0) 
                            continue;

                        positions.push_back(base_data[0]);       // ID
                        positions.push_back(base_data[1] + dx);  // Shifted X
                        positions.push_back(base_data[2] + dy);  // Shifted Y
                        positions.push_back(base_data[3] + dz);  // Shifted Z
                        positions.push_back(base_data[4]);       // Axis X
                        positions.push_back(base_data[5]);       // Axis Y
                        positions.push_back(base_data[6]);       // Axis Z
                        positions.push_back(base_data[7]);       // Angle
                    }
                }
            }
        }
    }

    return positions;
}

std::vector<double> PackingGeometry::sphere_array() const
{
    const std::size_t n = num_spheres();
    std::vector<double> sphere_positions;
    sphere_positions.reserve(n * 4);
    
    const auto& [sx, sy, sz, sr] = spheres;
    
    for (std::size_t i = 0; i < n; ++i)
    {
        sphere_positions.push_back(sx[i]);
        sphere_positions.push_back(sy[i]);
        sphere_positions.push_back(sz[i]);
        sphere_positions.push_back(sr[i]);
    }
    
    return sphere_positions;
}

void PackingGeometry::clear()
{
    // Clear all particle data arrays
    for (auto& vec : particles)
    {
        vec.clear();
    }
    
    // Clear all sphere data arrays
    for (auto& vec : spheres)
    {
        vec.clear();
    }
    
    // Reset particle offsets to initial state
    particle_offsets.clear();
    particle_offsets.push_back(0);
}

} // namespace pygrain3d

