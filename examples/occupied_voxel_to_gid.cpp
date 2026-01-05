#include "biomesh/VoxelMeshGenerator.hpp"
#include "biomesh/GiDExporter.hpp"
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomBuilder.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using namespace biomesh;

void printUsage(const char* programName) {
    std::cout << "\nUsage: " << programName << " <pdb_file> <voxel_size> <output_file> [padding]\n\n";
    std::cout << "Arguments:\n";
    std::cout << "  pdb_file     : Path to input PDB file\n";
    std::cout << "  voxel_size   : Edge length of voxels in Angstroms (e.g., 1.0)\n";
    std::cout << "  output_file  : Output GiD mesh file path (e.g., occupied_mesh.msh)\n";
    std::cout << "  padding      : Optional padding around bounding box in Angstroms (default: 2.0)\n\n";
    std::cout << "Example:\n";
    std::cout << "  " << programName << " protein.pdb 1.0 occupied_mesh.msh 2.0\n\n";
    std::cout << "Note: Occupied voxel meshes represent the molecule volume.\n";
    std::cout << "      Use larger voxel sizes for initial testing.\n\n";
}

int main(int argc, char* argv[]) {
    std::cout << "BioMesh - Occupied Voxel Mesh Generator with GiD Export\n";
    std::cout << "========================================================\n\n";
    
    // Check arguments
    if (argc < 4) {
        std::cerr << "Error: Insufficient arguments\n";
        printUsage(argv[0]);
        return 1;
    }
    
    std::string pdbFile = argv[1];
    double voxelSize = 0.0;
    std::string outputFile = argv[3];
    double padding = 2.0; // Default padding
    
    // Parse voxel size
    try {
        voxelSize = std::stod(argv[2]);
        if (voxelSize <= 0.0) {
            std::cerr << "Error: Voxel size must be positive\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid voxel size: " << argv[2] << "\n";
        return 1;
    }
    
    // Parse optional padding
    if (argc > 4) {
        try {
            padding = std::stod(argv[4]);
            if (padding < 0.0) {
                std::cerr << "Error: Padding cannot be negative\n";
                return 1;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid padding value: " << argv[4] << "\n";
            return 1;
        }
    }
    
    try {
        // Step 1: Parse PDB file
        std::cout << "Loading PDB file: " << pdbFile << "\n";
        auto basicAtoms = PDBParser::parsePDBFile(pdbFile);
        std::cout << "  Loaded " << basicAtoms.size() << " atoms\n\n";
        
        // Step 2: Enrich atoms with physical properties
        std::cout << "Enriching atoms with physical properties...\n";
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(basicAtoms);
        std::cout << "  Enriched " << enrichedAtoms.size() << " atoms\n\n";
        
        // Step 3: Create voxel grid
        std::cout << "Creating voxel grid...\n";
        std::cout << "  Voxel size: " << voxelSize << " Å\n";
        std::cout << "  Padding: " << padding << " Å\n";
        VoxelGrid voxelGrid(enrichedAtoms, voxelSize, padding);
        
        std::cout << "\n";
        voxelGrid.printStatistics();
        std::cout << "\n";
        
        // Step 4: Generate hexahedral mesh from occupied voxels
        // Note: VoxelMeshGenerator::generateHexMesh() creates mesh elements only for
        // voxels that intersect with atom spheres (occupied voxels), representing
        // the molecule volume. For void space mesh, use EmptyVoxelMeshGenerator.
        std::cout << "Generating hexahedral mesh from occupied voxels...\n";
        HexMesh mesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);
        
        std::cout << "  Generated mesh:\n";
        std::cout << "    Nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "    Elements: " << mesh.getElementCount() << "\n";
        
        // Calculate mesh statistics
        int occupiedVoxelCount = voxelGrid.getOccupiedVoxelCount();
        int theoreticalNodes = occupiedVoxelCount * 8;
        double efficiency = 0.0;
        if (theoreticalNodes > 0) {
            efficiency = (1.0 - (double)mesh.getNodeCount() / theoreticalNodes) * 100.0;
        }
        std::cout << "    Node sharing efficiency: " << std::fixed << std::setprecision(1) 
                  << efficiency << "%\n\n";
        
        // Warn about large meshes
        if (mesh.getElementCount() > 100000) {
            std::cout << "WARNING: Large mesh detected (" << mesh.getElementCount() 
                      << " elements). File may be large.\n\n";
        }
        
        // Step 5: Export to GiD format
        std::cout << "Exporting to GiD format: " << outputFile << "\n";
        bool success = GiDExporter::exportToGiD(mesh, outputFile);
        
        if (success) {
            std::cout << "  Export successful!\n";
            std::cout << "\nMesh file written to: " << outputFile << "\n";
            std::cout << "You can now open this file in GiD or any compatible FEM/CFD software.\n";
            return 0;
        } else {
            std::cerr << "  Export failed!\n";
            return 1;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
