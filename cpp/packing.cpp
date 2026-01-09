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

    auto& [px, py, pz, pr, pid] = geometry_.particles;
    auto& [sx, sy, sz, sr] = geometry_.spheres;

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

} // namespace pygrain