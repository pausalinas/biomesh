#include "biomesh/GiDExporter.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace biomesh {

bool GiDExporter::exportToGiD(const HexMesh& mesh, const std::string& filename) {
    // Check if mesh is empty
    if (mesh.nodes.empty() || mesh.elements.empty()) {
        std::cerr << "Error: Cannot export empty mesh to GiD format\n";
        return false;
    }
    
    // Open output file
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << filename << "\n";
        return false;
    }
    
    // Write GiD header
    outFile << "MESH dimension 3 ElemType Hexahedra Nnode 8\n\n";
    
    // Write coordinates section
    outFile << "Coordinates\n";
    outFile << std::fixed << std::setprecision(6);
    
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const Point3D& node = mesh.nodes[i];
        // GiD uses 1-based indexing
        outFile << (i + 1) << " " << node.x << " " << node.y << " " << node.z << "\n";
    }
    
    outFile << "End Coordinates\n\n";
    
    // Write elements section
    outFile << "Elements\n";
    
    for (size_t i = 0; i < mesh.elements.size(); ++i) {
        const auto& element = mesh.elements[i];
        // GiD uses 1-based indexing
        outFile << (i + 1);
        for (int nodeIdx : element) {
            outFile << " " << (nodeIdx + 1);
        }
        outFile << "\n";
    }
    
    outFile << "End Elements\n";
    
    outFile.close();
    
    return true;
}

} // namespace biomesh
