#include "biomesh/BoundingBox.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace biomesh {

BoundingBox::BoundingBox(const std::vector<std::unique_ptr<Atom>>& atoms, double padding) {
    if (atoms.empty()) {
        throw std::runtime_error("Cannot create bounding box from empty atom vector");
    }
    
    calculateBounds(atoms, padding);
}

void BoundingBox::calculateBounds(const std::vector<std::unique_ptr<Atom>>& atoms, double padding) {
    // Initialize with extreme values
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double minZ = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    double maxZ = std::numeric_limits<double>::lowest();

    // Find bounds considering van der Waals radii
    for (const auto& atom : atoms) {
        double x = atom->getX();
        double y = atom->getY();
        double z = atom->getZ();
        double radius = atom->getAtomicRadius();

        // Expand coordinates by van der Waals radius to get sphere boundaries
        minX = std::min(minX, x - radius);
        minY = std::min(minY, y - radius);
        minZ = std::min(minZ, z - radius);
        maxX = std::max(maxX, x + radius);
        maxY = std::max(maxY, y + radius);
        maxZ = std::max(maxZ, z + radius);
    }

    // Apply additional padding
    min_ = Point3D(minX - padding, minY - padding, minZ - padding);
    max_ = Point3D(maxX + padding, maxY + padding, maxZ + padding);
}

Point3D BoundingBox::getCenter() const {
    return Point3D(
        (min_.x + max_.x) * 0.5,
        (min_.y + max_.y) * 0.5,
        (min_.z + max_.z) * 0.5
    );
}

Point3D BoundingBox::getDimensions() const {
    return Point3D(
        max_.x - min_.x,
        max_.y - min_.y,
        max_.z - min_.z
    );
}

double BoundingBox::getVolume() const {
    Point3D dimensions = getDimensions();
    return dimensions.x * dimensions.y * dimensions.z;
}

double BoundingBox::getSurfaceArea() const {
    Point3D dimensions = getDimensions();
    return 2.0 * (dimensions.x * dimensions.y + 
                  dimensions.y * dimensions.z + 
                  dimensions.z * dimensions.x);
}

std::vector<Point3D> BoundingBox::getCorners() const {
    return {
        Point3D(min_.x, min_.y, min_.z),  // 0: min corner
        Point3D(max_.x, min_.y, min_.z),  // 1: +x
        Point3D(min_.x, max_.y, min_.z),  // 2: +y
        Point3D(max_.x, max_.y, min_.z),  // 3: +x+y
        Point3D(min_.x, min_.y, max_.z),  // 4: +z
        Point3D(max_.x, min_.y, max_.z),  // 5: +x+z
        Point3D(min_.x, max_.y, max_.z),  // 6: +y+z
        Point3D(max_.x, max_.y, max_.z)   // 7: max corner
    };
}

bool BoundingBox::contains(const Point3D& point) const {
    return (point.x >= min_.x && point.x <= max_.x &&
            point.y >= min_.y && point.y <= max_.y &&
            point.z >= min_.z && point.z <= max_.z);
}

void BoundingBox::expandBy(double additionalPadding) {
    min_.x -= additionalPadding;
    min_.y -= additionalPadding;
    min_.z -= additionalPadding;
    max_.x += additionalPadding;
    max_.y += additionalPadding;
    max_.z += additionalPadding;
}

} // namespace biomesh