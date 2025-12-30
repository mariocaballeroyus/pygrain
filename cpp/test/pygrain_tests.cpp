#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "pygrain.hpp"
#include "particle_operations.hpp"

 /**
  * @brief Test adding a spherical particle to the packing.
  * 
  * Verifies that the sphere is added at the origin with correct radius.
  */
TEST_CASE("Add Sphere Particle", "[Packing]")
{
    constexpr double tol = 1e-12;

    pygrain::Packing packing;
    packing.add_sphere_particle(1.0);

    REQUIRE(packing.get_geometry().num_spheres() == 1);

    const auto& [x, y, z, radii] = packing.get_geometry().sphere_data;

    REQUIRE(x[0] == Approx(0.0).margin(tol));
    REQUIRE(y[0] == Approx(0.0).margin(tol));
    REQUIRE(z[0] == Approx(0.0).margin(tol));
    REQUIRE(radii[0] == Approx(1.0).margin(tol));
}

/**
 * @brief Test translating a particle.
 * 
 * Verifies that the sphere's position is updated correctly after translation.
 */
TEST_CASE("Translate Particle", "[ParticleGeometry]")
{
    constexpr double tol = 1e-12;

    pygrain::Packing packing;
    packing.add_sphere_particle(1.0);

    pygrain::translate_particle(packing.get_geometry(), 0, {1.0, 2.0, 3.0});

    const auto& [x, y, z, r] = packing.get_geometry().sphere_data;

    REQUIRE(x[0] == Approx(1.0).margin(tol));
    REQUIRE(y[0] == Approx(2.0).margin(tol));
    REQUIRE(z[0] == Approx(3.0).margin(tol));
}

/**
 * @brief Test rotating a particle.
 * 
 * Verifies that the sphere's position is updated correctly after rotation.
 */
TEST_CASE("Rotate Particle", "[ParticleGeometry]")
{
    constexpr double tol = 1e-12;

    pygrain::Packing packing;
    packing.add_sphere_particle(1.0);

    // Manually set the center to (1,0,0) for rotation
    auto& [px, py, pz, pr] = packing.get_geometry().particle_data;
    px[0] = 1.0;
    py[0] = 0.0;
    pz[0] = 0.0;

    const auto& [sx, sy, sz, sr] = packing.get_geometry().sphere_data;

    // 90 deg around Z-axis
    pygrain::rotate_particle(packing.get_geometry(), 0, {0.0, 0.0, -1.0}, M_PI/2);

    REQUIRE(sx[0] == Approx(1.0).margin(tol));
    REQUIRE(sy[0] == Approx(1.0).margin(tol));
    REQUIRE(sz[0] == Approx(0.0).margin(tol));

    // 90 deg around X-axis
    pygrain::rotate_particle(packing.get_geometry(), 0, {-1.0, 0.0, 0.0}, M_PI/2);
    REQUIRE(sx[0] == Approx(1.0).margin(tol));
    REQUIRE(sy[0] == Approx(0.0).margin(tol));
    REQUIRE(sz[0] == Approx(-1.0).margin(tol));

    // 270 deg around Y-axis
    pygrain::rotate_particle(packing.get_geometry(), 0, {0.0, 1.0, 0.0}, 3*M_PI/2);
    REQUIRE(sx[0] == Approx(2.0).margin(tol));
    REQUIRE(sy[0] == Approx(0.0).margin(tol));
    REQUIRE(sz[0] == Approx(0.0).margin(tol));

    // 180 deg around Z-axis
    pygrain::rotate_particle(packing.get_geometry(), 0, {0.0, 0.0, 1.0}, M_PI);
    REQUIRE(sx[0] == Approx(0.0).margin(tol));
    REQUIRE(sy[0] == Approx(0.0).margin(tol));
    REQUIRE(sz[0] == Approx(0.0).margin(tol));
}

/**
 * @brief Test sphere approximation of a spheroid.
 * 
 * Verifies that the spheres approximating the spheroid are internally tangent
 * to the spheroid surface.
 */
TEST_CASE("Sphere approximation of spheroid", "[Spheroid]")
{
    constexpr double tol = 1e-12;
    
    const double aspect_ratio = 2.0;
    const double minor_axis = 1.0;

    const double a = minor_axis / 2.0;
    const double c = a * aspect_ratio;

    // Create geometry and generate spheroid particle
    pygrain::PackingGeometry geometry;
    pygrain::generate_spheroid_particle(geometry, aspect_ratio, minor_axis);

    const auto& [sx, sy, sz, sr] = geometry.sphere_data;
    const std::size_t num_spheres = geometry.num_spheres();

    REQUIRE(num_spheres > 0);

    // Check that each sphere is internally tangent to the spheroid surface
    for (std::size_t i = 0; i < num_spheres; ++i)
    {
        const double R = sr[i];
        const double x0 = sx[i];

        // For a prolate spheroid x^2/c^2 + y^2/a^2 + z^2/a^2 = 1,
        // a sphere centered at (x0, 0, 0) with radius R is internally tangent
        // when the discriminant of the intersection equation is zero.
        const double A = 1.0 - (a * a) / (c * c);
        const double B = 2.0 * x0;
        const double C = x0 * x0 + a * a - R * R;

        const double discriminant = B * B - 4.0 * A * C;
        REQUIRE(discriminant == Approx(0.0).margin(tol));
    }
}


/**
 * @brief Test volume coverage of spheroid by approximating spheres.
 * 
 * Verifies that the spheres cover at least a specified fraction of the spheroid 
 * volume for various aspect ratios.
 */
TEST_CASE("Spheroid volume coverage", "[Spheroid]")
{
    constexpr double min_volume_coverage = 0.97;
    constexpr int N = 100;  // number of grid points per axis
    
    const double eq_diameter = 1.0;
    const double aspect_ratio = GENERATE(1.5, 2.0, 3.0, 5.0, 10.0);
    CAPTURE(aspect_ratio);

    // Compute spheroid semi-axes
    const double c = (eq_diameter / 2.0) * std::sqrt(aspect_ratio);
    const double a = (eq_diameter / 2.0) / std::sqrt(aspect_ratio);

    // Create geometry and generate spheroid particle
    pygrain::PackingGeometry geometry;
    pygrain::generate_spheroid_particle(geometry, aspect_ratio, eq_diameter);

    const auto& [sx, sy, sz, sr] = geometry.sphere_data;
    const std::size_t num_spheres = geometry.num_spheres();

    int pts_inside = 0;   // points inside spheroid
    int pts_covered = 0;  // points covered by spheres

    // Generate uniform mesh of points within the bounding box
    const double dx = 2.0 * c / (N - 1);
    const double dy = 2.0 * a / (N - 1);
    const double dz = 2.0 * a / (N - 1);

    for (int i = 0; i < N; ++i)
    {
        const double x = -c + i * dx;
        for (int j = 0; j < N; ++j)
        {
            const double y = -a + j * dy;
            for (int k = 0; k < N; ++k)
            {
                const double z = -a + k * dz;

                // Check if the point is inside the spheroid using the spheroid equation
                const bool inside_spheroid = (x * x / (c * c)) + 
                                             (y * y / (a * a)) + 
                                             (z * z / (a * a)) <= 1.0;

                if (!inside_spheroid)
                {
                    continue;
                }

                ++pts_inside;

                // Check if the point is inside any of the approximating spheres
                for (std::size_t idx = 0; idx < num_spheres; ++idx)
                {
                    const double R = sr[idx];
                    const double dist_squared = (x - sx[idx]) * (x - sx[idx]) +
                                                (y - sy[idx]) * (y - sy[idx]) +
                                                (z - sz[idx]) * (z - sz[idx]);
                    if (dist_squared <= R * R)
                    {
                        ++pts_covered;
                        break; // point is covered, no need to check other spheres
                    }
                }
            }
        }
    }

    const double volume_coverage = static_cast<double>(pts_covered) / 
                                   static_cast<double>(pts_inside);
    REQUIRE(volume_coverage >= min_volume_coverage);
}

/**
 * @brief Test volume coverage of cylinder by approximating spheres.
 * 
 * Verifies that the spheres cover at least a specified fraction of the cylinder 
 * volume for various aspect ratios.
 */
TEST_CASE("Cylinder volume coverage", "[Cylinder]")
{
    constexpr double min_volume_coverage = 0.84;
    constexpr int N = 100;  // number of grid points per axis
    
    const double aspect_ratio = GENERATE(1.0, 1.5, 2.0, 3.0, 5.0, 10.0);
    CAPTURE(aspect_ratio);
    const double diameter = 1.0;

    // Compute cylinder dimensions
    const double r = diameter / 2.0;
    const double L = aspect_ratio * diameter;
    const double half_L = L / 2.0;

    // Create geometry and generate cylinder particle
    pygrain::PackingGeometry geometry;
    pygrain::generate_cylinder_particle(geometry, aspect_ratio, diameter);

    const auto& [sx, sy, sz, sr] = geometry.sphere_data;
    const std::size_t num_spheres = geometry.num_spheres();

    int pts_inside = 0;   // points inside cylinder
    int pts_covered = 0;  // points covered by spheres

    // Generate uniform mesh of points within the bounding box
    const double dx = 2.0 * half_L / (N - 1);
    const double dy = 2.0 * r / (N - 1);
    const double dz = 2.0 * r / (N - 1);

    for (int i = 0; i < N; ++i)
    {
        const double x = -half_L + i * dx;
        for (int j = 0; j < N; ++j)
        {
            const double y = -r + j * dy;
            for (int k = 0; k < N; ++k)
            {
                const double z = -r + k * dz;

                // Check if the point is inside the cylinder
                const bool inside_cylinder = (std::abs(x) <= half_L) && (y * y + z * z <= r * r);

                if (!inside_cylinder)
                {
                    continue;
                }

                ++pts_inside;

                // Check if the point is inside any of the approximating spheres
                for (std::size_t s = 0; s < num_spheres; ++s)
                {
                    const double R = sr[s];
                    const double dist_squared = (x - sx[s]) * (x - sx[s]) +
                                                (y - sy[s]) * (y - sy[s]) +
                                                (z - sz[s]) * (z - sz[s]);
                    if (dist_squared <= R * R)
                    {
                        ++pts_covered;
                        break;  // Point is covered, no need to check other spheres
                    }
                }
            }
        }
    }

    const double volume_coverage = static_cast<double>(pts_covered) / 
                                   static_cast<double>(pts_inside);
    REQUIRE(volume_coverage >= min_volume_coverage);
}