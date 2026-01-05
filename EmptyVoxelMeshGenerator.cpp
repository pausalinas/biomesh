#include "biomesh/EmptyVoxelMeshGenerator.hpp"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace biomesh {

HexMesh EmptyVoxelMeshGenerator::generateHexMesh(const VoxelGrid& voxelGrid) {
    HexMesh mesh;
    
    const auto& emptyVoxels = voxelGrid.getEmptyVoxels();
    
    if (emptyVoxels.empty()) {
        return mesh; // Return empty mesh if no empty voxels
    }
    
    // Step 1: Compute corner nodes for all empty voxels
    // Pre-allocate vector for thread-safe indexed assignment
    std::vector<std::array<Point3D, 8>> allCornerNodes;
    allCornerNodes.resize(emptyVoxels.size());
    
    // Parallelize corner node computation with OpenMP
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (size_t i = 0; i < emptyVoxels.size(); ++i) {
        allCornerNodes[i] = computeCornerNodes(emptyVoxels[i]);
    }
    
    // Step 2: Assign unique node indices with deduplication
    std::unordered_map<Point3D, int, Point3DHash, Point3DEqual> nodeMap;
    mesh.elements = assignUniqueNodeIndices(allCornerNodes, nodeMap, mesh.nodes);
    
    return mesh;
}

std::array<Point3D, 8> EmptyVoxelMeshGenerator::computeCornerNodes(const Voxel& voxel) {
    // Standard hexahedral element node ordering:
    // 0: (min.x, min.y, min.z) - bottom-left-back
    // 1: (max.x, min.y, min.z) - bottom-right-back  
    // 2: (max.x, max.y, min.z) - bottom-right-front
    // 3: (min.x, max.y, min.z) - bottom-left-front
    // 4: (min.x, min.y, max.z) - top-left-back
    // 5: (max.x, min.y, max.z) - top-right-back
    // 6: (max.x, max.y, max.z) - top-right-front
    // 7: (min.x, max.y, max.z) - top-left-front
    
    return {{
        Point3D(voxel.min.x, voxel.min.y, voxel.min.z), // 0
        Point3D(voxel.max.x, voxel.min.y, voxel.min.z), // 1
        Point3D(voxel.max.x, voxel.max.y, voxel.min.z), // 2
        Point3D(voxel.min.x, voxel.max.y, voxel.min.z), // 3
        Point3D(voxel.min.x, voxel.min.y, voxel.max.z), // 4
        Point3D(voxel.max.x, voxel.min.y, voxel.max.z), // 5
        Point3D(voxel.max.x, voxel.max.y, voxel.max.z), // 6
        Point3D(voxel.min.x, voxel.max.y, voxel.max.z)  // 7
    }};
}

std::vector<std::array<int, 8>> EmptyVoxelMeshGenerator::assignUniqueNodeIndices(
    const std::vector<std::array<Point3D, 8>>& cornerNodes,
    std::unordered_map<Point3D, int, Point3DHash, Point3DEqual>& nodeMap,
    std::vector<Point3D>& uniqueNodes) {
    
    std::vector<std::array<int, 8>> elements;
    elements.reserve(cornerNodes.size());
    
    int nextNodeIndex = 0;
    
    for (const auto& elementCorners : cornerNodes) {
        std::array<int, 8> elementConnectivity;
        
        // For each corner of this element
        for (size_t i = 0; i < 8; ++i) {
            const Point3D& corner = elementCorners[i];
            
            // Check if this node already exists
            auto it = nodeMap.find(corner);
            if (it != nodeMap.end()) {
                // Node already exists, reuse its index
                elementConnectivity[i] = it->second;
            } else {
                // New node, assign it a new index
                elementConnectivity[i] = nextNodeIndex;
                nodeMap[corner] = nextNodeIndex;
                uniqueNodes.push_back(corner);
                ++nextNodeIndex;
            }
        }
        
        elements.push_back(elementConnectivity);
    }
    
    return elements;
}

} // namespace biomesh
