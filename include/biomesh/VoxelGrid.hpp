#pragma once

#include "biomesh/BoundingBox.hpp"
#include "biomesh/Atom.hpp"
#include <vector>
#include <memory>
#include <array>

namespace biomesh {

/**
 * @brief Index structure for 3D voxel grid
 */
struct VoxelIndex {
    int i, j, k;  // Indices in x, y, z directions
    
    VoxelIndex() : i(0), j(0), k(0) {}
    VoxelIndex(int i_val, int j_val, int k_val) : i(i_val), j(j_val), k(k_val) {}
    
    bool operator==(const VoxelIndex& other) const {
        return i == other.i && j == other.j && k == other.k;
    }
};

/**
 * @brief Single voxel in the grid
 */
struct Voxel {
    VoxelIndex index;           // Grid indices
    Point3D center;             // Center point of the voxel
    Point3D min;                // Minimum corner
    Point3D max;                // Maximum corner
    bool occupied;              // True if voxel intersects with any atom
    std::vector<size_t> atomIds; // IDs of atoms that intersect this voxel
    
    Voxel() : occupied(false) {}
};

/**
 * @brief Uniform voxel grid for molecular structures
 * 
 * Tessellates a bounding box into uniform cubic voxels and classifies
 * each voxel as occupied or empty based on intersection with atom spheres.
 */
class VoxelGrid {
public:
    /**
     * @brief Constructor with atoms and voxel size
     * @param atoms Vector of fully initialized atoms
     * @param voxelSize Edge length of each cubic voxel
     * @param padding Additional padding to add to the bounding box (default: 0.0)
     */
    VoxelGrid(const std::vector<std::unique_ptr<Atom>>& atoms, 
              double voxelSize,
              double padding = 0.0);
    
    /**
     * @brief Constructor with existing bounding box and voxel size
     * @param boundingBox Pre-computed bounding box
     * @param atoms Vector of fully initialized atoms
     * @param voxelSize Edge length of each cubic voxel
     */
    VoxelGrid(const BoundingBox& boundingBox,
              const std::vector<std::unique_ptr<Atom>>& atoms,
              double voxelSize);
    
    /**
     * @brief Get voxel size (edge length)
     */
    double getVoxelSize() const { return voxelSize_; }
    
    /**
     * @brief Get grid dimensions (number of voxels in each direction)
     */
    const std::array<int, 3>& getDimensions() const { return dimensions_; }
    
    /**
     * @brief Get total number of voxels
     */
    int getTotalVoxelCount() const { return totalVoxels_; }
    
    /**
     * @brief Get number of occupied voxels
     */
    int getOccupiedVoxelCount() const { return occupiedVoxels_.size(); }
    
    /**
     * @brief Get number of empty voxels
     */
    int getEmptyVoxelCount() const { return emptyVoxels_.size(); }
    
    /**
     * @brief Get all occupied voxels
     */
    const std::vector<Voxel>& getOccupiedVoxels() const { return occupiedVoxels_; }
    
    /**
     * @brief Get all empty voxels
     */
    const std::vector<Voxel>& getEmptyVoxels() const { return emptyVoxels_; }
    
    /**
     * @brief Get voxel at specific grid indices
     * @param i Index in x direction
     * @param j Index in y direction
     * @param k Index in z direction
     * @return Pointer to voxel, or nullptr if indices are out of bounds
     */
    const Voxel* getVoxel(int i, int j, int k) const;
    
    /**
     * @brief Get bounding box used for voxelization
     */
    const BoundingBox& getBoundingBox() const { return boundingBox_; }
    
    /**
     * @brief Print grid statistics
     */
    void printStatistics() const;

private:
    /**
     * @brief Initialize the voxel grid
     * @param atoms Vector of atoms for occupancy detection
     */
    void initialize(const std::vector<std::unique_ptr<Atom>>& atoms);
    
    /**
     * @brief Check if a voxel center intersects with any atom sphere
     * @param voxelCenter Center point of the voxel
     * @param atoms Vector of atoms to check
     * @param intersectingAtomIds Output vector for atom IDs that intersect
     * @return true if voxel is occupied
     */
    bool isVoxelOccupied(const Point3D& voxelCenter,
                        const std::vector<std::unique_ptr<Atom>>& atoms,
                        std::vector<size_t>& intersectingAtomIds) const;
    
    /**
     * @brief Convert grid indices to linear index
     * @param i Index in x direction
     * @param j Index in y direction
     * @param k Index in z direction
     * @return Linear index
     */
    int getLinearIndex(int i, int j, int k) const {
        return i + j * dimensions_[0] + k * dimensions_[0] * dimensions_[1];
    }
    
    BoundingBox boundingBox_;              // Bounding box for the grid
    double voxelSize_;                     // Edge length of each voxel
    std::array<int, 3> dimensions_;        // Grid dimensions (nx, ny, nz)
    int totalVoxels_;                      // Total number of voxels
    
    std::vector<Voxel> occupiedVoxels_;    // Voxels that intersect atoms
    std::vector<Voxel> emptyVoxels_;       // Voxels that don't intersect atoms
    std::vector<Voxel*> voxelGrid_;        // 3D grid for fast lookup (nullptr for empty)
};

} // namespace biomesh
