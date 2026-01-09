#ifndef PYPACKING_HPP
#define PYPACKING_HPP

#include "geometry.hpp"

enum class BoundaryType
{
    PERIODIC,
    WALL
};

class Packing
{
public:
    explicit Packing(std::array<double, 3> lengths_)
    : boundary_types_({BoundaryType::PERIODIC, 
                       BoundaryType::PERIODIC, 
                       BoundaryType::PERIODIC}),
      lengths_(lengths_) {}

    Packing(std::array<double, 3> lengths,
            std::array<BoundaryType, 3> boundary_types)
    : boundary_types_(boundary_types), lengths_(lengths) {}
    

private:
    std::array<BoundaryType, 3> boundary_types_;
    std::array<double, 3> lengths_;
    PackingGeometry geometry_;
};

#endif // PYPACKING_HPP