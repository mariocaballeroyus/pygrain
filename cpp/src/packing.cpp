#include "packing.hpp"
#include "particle_operations.hpp"

#include <algorithm>
#include <cmath>

namespace pygrain
{

void Packing::add_sphere_particle(double radius)
{
    generate_sphere_particle(geometry_, radius);
}

void Packing::add_spheroid_particle(double aspect_ratio, double eq_diameter)
{
    generate_spheroid_particle(geometry_, aspect_ratio, eq_diameter);
}

void Packing::add_cylinder_particle(double aspect_ratio, double eq_diameter)
{
    generate_cylinder_particle(geometry_, aspect_ratio, eq_diameter);
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

void Packing::generate(unsigned int max_iterations)
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

    auto& [px, py, pz, pr] = geometry_.particle_data;
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
        std::cout << "Iteration " << num_iterations 
                  << ": Number of overlaps = " << num_overlaps << std::endl;
    }
    while (num_overlaps > 0 && num_iterations < max_iterations);
}

void Packing::export_spheres_csv(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    file << "x,y,z,r,particle_id" << std::endl;

    const auto& [px, py, pz, pr] = geometry_.particle_data;
    const auto& [sx, sy, sz, sr] = geometry_.sphere_data;
    
    const auto& offsets = geometry_.particle_offsets;

    std::size_t num_exported = 0;

    for (std::size_t p = 0; p < geometry_.num_particles(); ++p)
    {
        const double p_cx = px[p];
        const double p_cy = py[p];
        const double p_cz = pz[p];
        const double p_bound = pr[p];

        const std::size_t start = offsets[p];
        const std::size_t end = offsets[p + 1];

        // Check all 27 potential image locations (including the original box dx=dy=dz=0)
        for (int dx = -1; dx <= 1; ++dx)
        {
            for (int dy = -1; dy <= 1; ++dy)
            {
                for (int dz = -1; dz <= 1; ++dz)
                {
                    // Center of the particle in this specific ghost image
                    const double tx_p = p_cx + dx * lengths_[0];
                    const double ty_p = p_cy + dy * lengths_[1];
                    const double tz_p = p_cz + dz * lengths_[2];

                    // Check if the particle's bounding sphere touches the primary box [0, L]
                    // This uses a simple AABB vs Sphere intersection logic
                    if (tx_p + p_bound > 0 && tx_p - p_bound < lengths_[0] &&
                        ty_p + p_bound > 0 && ty_p - p_bound < lengths_[1] &&
                        tz_p + p_bound > 0 && tz_p - p_bound < lengths_[2])
                    {
                        // The particle (or at least its bounding sphere) is visible here.
                        // Export all constituent spheres for this image.
                        for (std::size_t s = start; s < end; ++s)
                        {
                            file << sx[s] + dx * lengths_[0] << ","
                                 << sy[s] + dy * lengths_[1] << ","
                                 << sz[s] + dz * lengths_[2] << ","
                                 << sr[s] << ","
                                 << p << std::endl;
                            ++num_exported;
                        }
                    }
                }
            }
        }
    }

    file.close();
    std::cout << "Exported " << num_exported << " spheres (including periodic images) to " << filename << std::endl;
}

} // namespace pygrain