#include "pygrain.hpp"

int main()
{
    // Create packing with specified domain size
    std::array<double, 3> domain = {100.0, 100.0, 100.0};
    pygrain::Packing packing(domain);
    
    // Add cylindrical particles
    const int num_particles = 100;
    const double aspect_ratio = 2.0;   // length / diameter
    const double diameter = 15.0;
    const double length = aspect_ratio * diameter;
    
    // Cylinder volume: V = pi * r^2 * L
    const double cylinder_volume = M_PI * (diameter / 2.0) * (diameter / 2.0) * length;
    double solid_volume = 0.0;
    
    // Add particles and accumulate solid volume
    for (int i = 0; i < num_particles; ++i)
    {
        packing.add_cylinder_particle(aspect_ratio, diameter);
        solid_volume += cylinder_volume;
    }
    
    // Compute porosity
    const double domain_volume = domain[0] * domain[1] * domain[2];
    const double porosity = 1.0 - (solid_volume / domain_volume);
    
    std::cout << "Created " << packing.get_geometry().num_particles() << " cylinder particles" << std::endl;
    std::cout << "Total spheres: " << packing.get_geometry().num_spheres() << std::endl;
    std::cout << "Porosity: " << porosity << std::endl;
    
    // Generate packing (resolve overlaps)
    const unsigned int max_iterations = 10000;
    packing.generate(max_iterations);
    std::cout << "Generated packing with porosity: " << porosity << std::endl;
    
    // Export to CSV
    packing.export_spheres_csv("packing.csv");
    
    return 0;
}