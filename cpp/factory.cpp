#include "geometry.hpp"

#include <cmath>
#include <algorithm>

namespace pygrain
{

void begin_particle(PackingGeometry& geometry, 
                    double x, double y, double z, 
                    double bounding_radius,
                    std::size_t geometry_idx)
{
    auto& [px, py, pz, pr, pid] = geometry.particle_data;

    px.push_back(x);
    py.push_back(y);
    pz.push_back(z);
    pr.push_back(bounding_radius);
    pid.push_back(geometry_idx);
}

void end_particle(PackingGeometry& geometry)
{
    geometry.particle_offsets.push_back(geometry.sphere_data[0].size());
}

void add_particle_sphere(PackingGeometry& geometry, 
                         double x, double y, double z, 
                         double radius)
{
    auto& [sx, sy, sz, sr] = geometry.sphere_data;

    sx.push_back(x);
    sy.push_back(y);
    sz.push_back(z);
    sr.push_back(radius);
}

void generate_sphere_particle(PackingGeometry& geometry, 
                              double radius,
                              std::size_t geometry_idx)
{
    begin_particle(geometry, 0.0, 0.0, 0.0, radius, geometry_idx);
    add_particle_sphere(geometry, 0.0, 0.0, 0.0, radius);
    end_particle(geometry);
}

void generate_spheroid_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double minor_axis,
                                std::size_t geometry_idx, 
                                double sphere_precision)
{
    const double a = minor_axis / 2.0;  // minor semi-axis
    const double c = a * aspect_ratio;       // major semi-axis

    begin_particle(geometry, 0.0, 0.0, 0.0, c, geometry_idx);

    // If nearly spherical, just add one sphere
    if (std::abs(a - c) / std::max(a, c) < 1e-4)
    {
        add_particle_sphere(geometry, 0.0, 0.0, 0.0, a);
        end_particle(geometry);
        return;
    }

    double x0, x, y, radius;
    double last_radius = 0.0;

    x = c;
    while (true)
    {
        x0 = x * (1.0 - (a * a) / (c * c));
        y = a * std::sqrt(1.0 - (x * x) / (c * c));
        radius = std::sqrt((x - x0) * (x - x0) + y * y);

        // Add spheres symmetrically along the major axis (z-axis for prolate)
        add_particle_sphere(geometry, x0, 0.0, 0.0, radius);
        if (x0 > 1e-10)  // Only add second sphere if not at center
        {
            add_particle_sphere(geometry, -x0, 0.0, 0.0, radius);
        }

        last_radius = radius;
        x -= radius / sphere_precision;

        if (x <= 0.0)
            break;
    }

    // Add center sphere if the last radius is smaller than minor axis
    if (last_radius <= 0.97 * a)
    {
        add_particle_sphere(geometry, 0.0, 0.0, 0.0, a);
    }

    end_particle(geometry);
}

void generate_cylinder_particle(PackingGeometry& geometry, 
                                double aspect_ratio, 
                                double diameter,
                                std::size_t geometry_idx, 
                                double sphere_precision)
{    
    const double r = diameter / 2.0;           // cylinder radius
    const double L = aspect_ratio * diameter;  // cylinder length
    const double half_L = L / 2.0;

    const double bounding_radius = std::sqrt(half_L * half_L + r * r);
    begin_particle(geometry, 0.0, 0.0, 0.0, bounding_radius, geometry_idx);

    // Step size for placing spheres along axis
    const double step = r / sphere_precision;

    // Add large spheres (radius = r) along the central axis
    const double end_sphere_x = std::max(0.0, half_L - r);
    
    double x = 0.0;
    while (x <= end_sphere_x + 1e-10)
    {
        add_particle_sphere(geometry, x, 0.0, 0.0, r);
        if (x > 1e-10)
        {
            add_particle_sphere(geometry, -x, 0.0, 0.0, r);
        }
        x += step;
    }
    
    if (end_sphere_x > 0 && std::abs(x - step - end_sphere_x) > 1e-10)
    {
        add_particle_sphere(geometry, end_sphere_x, 0.0, 0.0, r);
        add_particle_sphere(geometry, -end_sphere_x, 0.0, 0.0, r);
    }

    // Add two rings of small spheres at the corners (cap edges)
    const double edge_radius = r / (2.0 * sphere_precision);
    const int num_edge_spheres = static_cast<int>(std::ceil(M_PI * r / edge_radius));
    const double edge_angle_step = 2.0 * M_PI / num_edge_spheres;

    // Outer ring: at the very edge where cap meets curved surface
    const double outer_radial = r - edge_radius;
    const double outer_x = half_L - edge_radius;

    for (int i = 0; i < num_edge_spheres; ++i)
    {
        const double angle = i * edge_angle_step;
        const double y = outer_radial * std::cos(angle);
        const double z = outer_radial * std::sin(angle);

        add_particle_sphere(geometry, outer_x, y, z, edge_radius);
        add_particle_sphere(geometry, -outer_x, y, z, edge_radius);
    }

    // Inner ring: at half the radial distance to fill the gap
    const double inner_radial = r * 0.5;
    const double inner_x = half_L - edge_radius;
    const int num_inner_spheres = static_cast<int>(std::ceil(M_PI * inner_radial / edge_radius));
    const double inner_angle_step = 2.0 * M_PI / num_inner_spheres;

    for (int i = 0; i < num_inner_spheres; ++i)
    {
        const double angle = i * inner_angle_step;
        const double y = inner_radial * std::cos(angle);
        const double z = inner_radial * std::sin(angle);

        add_particle_sphere(geometry, inner_x, y, z, edge_radius);
        add_particle_sphere(geometry, -inner_x, y, z, edge_radius);
    }

    end_particle(geometry);
}

} // namespace pygrain