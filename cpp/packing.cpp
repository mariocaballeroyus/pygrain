#include "packing.hpp"
#include "operations.hpp"

#include <algorithm>
#include <cmath>

namespace pygrain3d
{

void Packing::add_sphere_particles(double radius, 
                                   int num,
                                   std::size_t geometry_idx)
{
    for (int i = 0; i < num; ++i)
    {
        generate_sphere_particle(geometry_, radius, geometry_idx);
    }
}

void Packing::add_spheroid_particles(double aspect_ratio, 
                                     double eq_diameter, 
                                     int num,
                                     std::size_t geometry_idx)
{
    for (int i = 0; i < num; ++i)
    {
        generate_spheroid_particle(geometry_, aspect_ratio, eq_diameter, geometry_idx);
    }
}

void Packing::add_cylinder_particles(double aspect_ratio, 
                                     double eq_diameter, 
                                     int num,
                                     std::size_t geometry_idx)
{
    for (int i = 0; i < num; ++i)
    {
        generate_cylinder_particle(geometry_, aspect_ratio, eq_diameter, geometry_idx);
    }
}

void Packing::randomize_particles()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());

    std::uniform_real_distribution<double> dist_x(0.0, lengths_[0]),
                                           dist_y(0.0, lengths_[1]),
                                           dist_z(0.0, lengths_[2]),
                                           dist_angle(0.0, 2.0 * M_PI),
                                           dist_cos(-1.0, 1.0);

    for (std::size_t i = 0; i < geometry_.num_particles(); ++i)
    {
        std::array<double, 3> random_position = {dist_x(gen),
                                                 dist_y(gen),
                                                 dist_z(gen)};

        double theta = dist_angle(gen);   // azimuthal angle
        double cos_phi = dist_cos(gen);   // cos(phi) uniformly in [-1,1]
        double phi = std::acos(cos_phi);  // polar angle

        std::array<double, 3> random_axis = {std::sin(phi) * std::cos(theta),
                                             std::sin(phi) * std::sin(theta),
                                             std::cos(phi)};
        double random_angle = dist_angle(gen);

        set_particle_position(geometry_, i, random_position);
        rotate_particle(geometry_, i, random_axis, random_angle);
    }
}

void Packing::generate(unsigned int max_iterations,
                       unsigned int log_interval)
{
    randomize_particles();

    unsigned int num_iterations = 0;
    std::size_t num_overlaps;

    static std::random_device rd;
    static std::mt19937 gen(rd());

    // Vector of particle indices for shuffling
    std::vector<std::size_t> indices(geometry_.num_particles());
    for (std::size_t i = 0; i < indices.size(); ++i)
        indices[i] = i;

    auto& [px, py, pz, pr, pid] = geometry_.particle_data;
    auto& [sx, sy, sz, sr] = geometry_.sphere_data;

    do 
    {
        num_overlaps = 0;
        
        // Shuffle particle order to avoid bias
        std::shuffle(indices.begin(), indices.end(), gen);

        // Resolution pass
        for (std::size_t ii = 0; ii < indices.size(); ++ii)
        {
            const std::size_t i = indices[ii];

            for (std::size_t jj = ii + 1; jj < indices.size(); ++jj)
            {
                const std::size_t j = indices[jj];

                // Minimum image convention - compute shift to bring j closer to i
                std::array<double, 3> img_displacement = {0.0, 0.0, 0.0};
                const double delta_x = px[i] - px[j];
                const double delta_y = py[i] - py[j];
                const double delta_z = pz[i] - pz[j];

                if (delta_x > 0.5 * lengths_[0])
                    img_displacement[0] = lengths_[0];
                else if (delta_x < -0.5 * lengths_[0])
                    img_displacement[0] = -lengths_[0];

                if (delta_y > 0.5 * lengths_[1])
                    img_displacement[1] = lengths_[1];
                else if (delta_y < -0.5 * lengths_[1])
                    img_displacement[1] = -lengths_[1];

                if (delta_z > 0.5 * lengths_[2])
                    img_displacement[2] = lengths_[2];
                else if (delta_z < -0.5 * lengths_[2])
                    img_displacement[2] = -lengths_[2];

                // Temporarily shift j to its closest periodic image
                translate_particle(geometry_, j, img_displacement);

                const double dx_dist = px[i] - px[j];
                const double dy_dist = py[i] - py[j];
                const double dz_dist = pz[i] - pz[j];
                const double distance = std::sqrt(dx_dist*dx_dist + dy_dist*dy_dist + dz_dist*dz_dist);
                const double radius_sum = pr[i] + pr[j];

                if (distance < radius_sum) // potential overlap
                {
                    if (particle_interaction(geometry_, i, j))
                        ++num_overlaps;
                }

                // Shift j back (we only shifted temporarily for distance calculation)
                img_displacement[0] = -img_displacement[0];
                img_displacement[1] = -img_displacement[1];
                img_displacement[2] = -img_displacement[2];
                translate_particle(geometry_, j, img_displacement);
            }

            // Keep particle i inside the primary box
            std::array<double, 3> shift_inside = {0.0, 0.0, 0.0};
            if (px[i] < 0.0)
                shift_inside[0] = lengths_[0];
            else if (px[i] >= lengths_[0])
                shift_inside[0] = -lengths_[0];

            if (py[i] < 0.0)
                shift_inside[1] = lengths_[1];
            else if (py[i] >= lengths_[1])
                shift_inside[1] = -lengths_[1];

            if (pz[i] < 0.0)
                shift_inside[2] = lengths_[2];
            else if (pz[i] >= lengths_[2])
                shift_inside[2] = -lengths_[2];

            if (shift_inside[0] != 0.0 || shift_inside[1] != 0.0 || shift_inside[2] != 0.0)
                translate_particle(geometry_, i, shift_inside);
        }

        ++num_iterations;
        if (log_interval > 0 && num_iterations % log_interval == 0)
        {
        std::cout << "Iteration " << num_iterations 
                  << ": Number of overlaps = " << num_overlaps << std::endl;
        }
    }
    while (num_overlaps > 0 && num_iterations < max_iterations);

    std::cout << "Packing generation completed at iteration " << num_iterations << std::endl;
}

std::array<double, 8> Packing::particle_data(std::size_t idx) const
{
    const auto& [px, py, pz, pr, pid] = geometry_.particle_data;
    const auto& [sx, sy, sz, sr] = geometry_.sphere_data;
    const auto& offsets = geometry_.particle_offsets;

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

std::vector<double> Packing::data_array(bool periodic) const
{
    const std::size_t n = geometry_.num_particles();
    std::vector<double> positions;
    
    // Optimization: Reserve space for at least the base particles
    positions.reserve(n * 8);

    // Add all base particles regardless of periodic value ---
    for (std::size_t i = 0; i < n; ++i)
    {
        auto pos = particle_data(i); 
        for (std::size_t j = 0; j < 8; ++j)
        {
            positions.push_back(pos[j]);
        }
    }

    // --- Step 2 & 3: Handle Periodic Images ---
    if (periodic)
    {
        const auto& [px, py, pz, pr, pid] = geometry_.particle_data;

        for (std::size_t i = 0; i < n; ++i)
        {
            const double radius = pr[i];
            
            // Determine possible shifts for each dimension
            // A shift is needed if the center is within 'radius' of a boundary
            std::vector<double> dx_shifts = {0.0};
            if (px[i] < radius) dx_shifts.push_back(lengths_[0]);
            else if (px[i] > lengths_[0] - radius) dx_shifts.push_back(-lengths_[0]);

            std::vector<double> dy_shifts = {0.0};
            if (py[i] < radius) dy_shifts.push_back(lengths_[1]);
            else if (py[i] > lengths_[1] - radius) dy_shifts.push_back(-lengths_[1]);

            std::vector<double> dz_shifts = {0.0};
            if (pz[i] < radius) dz_shifts.push_back(lengths_[2]);
            else if (pz[i] > lengths_[2] - radius) dz_shifts.push_back(-lengths_[2]);

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

} // namespace pygrain