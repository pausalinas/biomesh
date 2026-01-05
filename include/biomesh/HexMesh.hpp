#pragma once

#include "biomesh/BoundingBox.hpp" // For Point3D
#include <vector>
#include <array>
#include <cmath>
#include <unordered_map>

namespace biomesh {

/**
 * @brief Hash function for Point3D to enable use in std::unordered_map
 */
struct Point3DHash {
    std::size_t operator()(const Point3D& p) const {
        // Use a simple hash combination
        std::size_t h1 = std::hash<double>{}(p.x);
        std::size_t h2 = std::hash<double>{}(p.y);
        std::size_t h3 = std::hash<double>{}(p.z);
        
        // Combine hashes using a simple method
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

/**
 * @brief Equality function for Point3D to enable use in std::unordered_map
 */
struct Point3DEqual {
    bool operator()(const Point3D& lhs, const Point3D& rhs) const {
        const double epsilon = 1e-12; // Small tolerance for floating point comparison
        return std::abs(lhs.x - rhs.x) < epsilon &&
               std::abs(lhs.y - rhs.y) < epsilon &&
               std::abs(lhs.z - rhs.z) < epsilon;
    }
};

/**
 * @brief Hexahedral mesh data structure
 */
struct HexMesh {
    std::vector<Point3D> nodes;                     // Unique node coordinates
    std::vector<std::array<int, 8>> elements;       // Element connectivity (8 node indices per hex)
    
    /**
     * @brief Get number of nodes in the mesh
     */
    size_t getNodeCount() const { return nodes.size(); }
    
    /**
     * @brief Get number of elements in the mesh
     */
    size_t getElementCount() const { return elements.size(); }
};

} // namespace biomesh
