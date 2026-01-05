#pragma once

#include "biomesh/Atom.hpp"
#include <vector>
#include <memory>

namespace biomesh {

/**
 * @brief 3D point structure
 */
struct Point3D {
    double x, y, z;
    
    Point3D() : x(0.0), y(0.0), z(0.0) {}
    Point3D(double x_val, double y_val, double z_val) : x(x_val), y(y_val), z(z_val) {}
};

/**
 * @brief Bounding box calculator for molecular structures
 * 
 * Calculates the 3D bounding box that encompasses all atoms including
 * their van der Waals radii, with optional padding.
 */
class BoundingBox {
public:
    /**
     * @brief Constructor with atoms and optional padding
     * @param atoms Vector of fully initialized atoms
     * @param padding Additional padding to add to the bounding box (default: 0.0)
     */
    explicit BoundingBox(const std::vector<std::unique_ptr<Atom>>& atoms, double padding = 0.0);

    /**
     * @brief Get minimum corner of the bounding box
     * @return Point3D representing the minimum corner
     */
    const Point3D& getMin() const { return min_; }

    /**
     * @brief Get maximum corner of the bounding box
     * @return Point3D representing the maximum corner
     */
    const Point3D& getMax() const { return max_; }

    /**
     * @brief Get center point of the bounding box
     * @return Point3D representing the center
     */
    Point3D getCenter() const;

    /**
     * @brief Get dimensions of the bounding box
     * @return Point3D with x, y, z dimensions
     */
    Point3D getDimensions() const;

    /**
     * @brief Get volume of the bounding box
     * @return Volume in cubic Angstroms
     */
    double getVolume() const;

    /**
     * @brief Get surface area of the bounding box
     * @return Surface area in square Angstroms
     */
    double getSurfaceArea() const;

    /**
     * @brief Get all 8 corner coordinates of the bounding box
     * @return Vector of 8 Point3D representing all corners
     */
    std::vector<Point3D> getCorners() const;

    /**
     * @brief Check if a point is inside the bounding box
     * @param point Point to check
     * @return true if point is inside the box
     */
    bool contains(const Point3D& point) const;

    /**
     * @brief Expand the bounding box by additional padding
     * @param additionalPadding Additional padding to add
     */
    void expandBy(double additionalPadding);

private:
    /**
     * @brief Calculate bounding box from atoms
     * @param atoms Vector of atoms
     * @param padding Padding to add
     */
    void calculateBounds(const std::vector<std::unique_ptr<Atom>>& atoms, double padding);

    Point3D min_;
    Point3D max_;
};

} // namespace biomesh