#include "biomesh/VoxelGrid.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace biomesh {

VoxelGrid::VoxelGrid(const std::vector<std::unique_ptr<Atom>>& atoms, 
                     double voxelSize,
                     double padding)
    : boundingBox_(atoms, padding), voxelSize_(voxelSize) {
    
    if (voxelSize <= 0.0) {
        throw std::invalid_argument("Voxel size must be positive");
    }
    
    initialize(atoms);
}

VoxelGrid::VoxelGrid(const BoundingBox& boundingBox,
                     const std::vector<std::unique_ptr<Atom>>& atoms,
                     double voxelSize)
    : boundingBox_(boundingBox), voxelSize_(voxelSize) {
    
    if (voxelSize <= 0.0) {
        throw std::invalid_argument("Voxel size must be positive");
    }
    
    initialize(atoms);
}

void VoxelGrid::initialize(const std::vector<std::unique_ptr<Atom>>& atoms) {
    // Calculate grid dimensions
    Point3D dims = boundingBox_.getDimensions();
    dimensions_[0] = static_cast<int>(std::ceil(dims.x / voxelSize_));
    dimensions_[1] = static_cast<int>(std::ceil(dims.y / voxelSize_));
    dimensions_[2] = static_cast<int>(std::ceil(dims.z / voxelSize_));
    
    // Ensure at least one voxel in each direction
    dimensions_[0] = std::max(1, dimensions_[0]);
    dimensions_[1] = std::max(1, dimensions_[1]);
    dimensions_[2] = std::max(1, dimensions_[2]);
    
    totalVoxels_ = dimensions_[0] * dimensions_[1] * dimensions_[2];
    
    // Initialize voxel grid for fast lookup
    voxelGrid_.resize(totalVoxels_, nullptr);
    
    // Reserve space for voxels (estimate 30% occupancy)
    occupiedVoxels_.reserve(totalVoxels_ / 3);
    emptyVoxels_.reserve(totalVoxels_);
    
    Point3D min = boundingBox_.getMin();
    
    // Generate all voxels and classify them
    for (int k = 0; k < dimensions_[2]; ++k) {
        for (int j = 0; j < dimensions_[1]; ++j) {
            for (int i = 0; i < dimensions_[0]; ++i) {
                Voxel voxel;
                voxel.index = VoxelIndex(i, j, k);
                
                // Calculate voxel bounds
                voxel.min.x = min.x + i * voxelSize_;
                voxel.min.y = min.y + j * voxelSize_;
                voxel.min.z = min.z + k * voxelSize_;
                
                voxel.max.x = voxel.min.x + voxelSize_;
                voxel.max.y = voxel.min.y + voxelSize_;
                voxel.max.z = voxel.min.z + voxelSize_;
                
                // Calculate voxel center
                voxel.center.x = voxel.min.x + voxelSize_ * 0.5;
                voxel.center.y = voxel.min.y + voxelSize_ * 0.5;
                voxel.center.z = voxel.min.z + voxelSize_ * 0.5;
                
                // Check occupancy
                voxel.occupied = isVoxelOccupied(voxel.center, atoms, voxel.atomIds);
                
                // Store in appropriate vector
                if (voxel.occupied) {
                    occupiedVoxels_.push_back(voxel);
                    voxelGrid_[getLinearIndex(i, j, k)] = &occupiedVoxels_.back();
                } else {
                    emptyVoxels_.push_back(voxel);
                    // Note: we could also store empty voxels in the grid if needed
                }
            }
        }
    }
}

bool VoxelGrid::isVoxelOccupied(const Point3D& voxelCenter,
                                const std::vector<std::unique_ptr<Atom>>& atoms,
                                std::vector<size_t>& intersectingAtomIds) const {
    intersectingAtomIds.clear();
    
    for (const auto& atom : atoms) {
        // Calculate distance from voxel center to atom center
        double dx = voxelCenter.x - atom->getX();
        double dy = voxelCenter.y - atom->getY();
        double dz = voxelCenter.z - atom->getZ();
        double distanceSquared = dx * dx + dy * dy + dz * dz;
        
        // Check if voxel center is within atom sphere
        double radius = atom->getAtomicRadius();
        if (distanceSquared <= radius * radius) {
            intersectingAtomIds.push_back(atom->getId());
        }
    }
    
    return !intersectingAtomIds.empty();
}

const Voxel* VoxelGrid::getVoxel(int i, int j, int k) const {
    if (i < 0 || i >= dimensions_[0] ||
        j < 0 || j >= dimensions_[1] ||
        k < 0 || k >= dimensions_[2]) {
        return nullptr;
    }
    
    return voxelGrid_[getLinearIndex(i, j, k)];
}

void VoxelGrid::printStatistics() const {
    std::cout << "\n=== Voxel Grid Statistics ===\n";
    std::cout << "Voxel size: " << voxelSize_ << " Å\n";
    std::cout << "Grid dimensions: " << dimensions_[0] << " x " 
              << dimensions_[1] << " x " << dimensions_[2] << "\n";
    std::cout << "Total voxels: " << totalVoxels_ << "\n";
    std::cout << "Occupied voxels: " << occupiedVoxels_.size() 
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * occupiedVoxels_.size() / totalVoxels_) << "%)\n";
    std::cout << "Empty voxels: " << emptyVoxels_.size()
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * emptyVoxels_.size() / totalVoxels_) << "%)\n";
    
    Point3D dims = boundingBox_.getDimensions();
    std::cout << "Bounding box dimensions: " 
              << std::fixed << std::setprecision(3)
              << dims.x << " x " << dims.y << " x " << dims.z << " Å\n";
    std::cout << "Bounding box volume: " 
              << std::fixed << std::setprecision(2)
              << boundingBox_.getVolume() << " Ų\n";
    std::cout << "=============================\n\n";
}

} // namespace biomesh
