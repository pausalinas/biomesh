#include "biomesh/VoxelGrid.hpp"
#include "biomesh/VoxelMeshGenerator.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomBuilder.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using namespace biomesh;

void demonstrateVoxelization() {
    std::cout << "=== Voxelization Demo ===\n\n";
    
    // Create test atoms
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    
    // Create a simple molecular structure
    basicAtoms.push_back(std::make_unique<Atom>("N"));
    basicAtoms.back()->setCoordinates(20.0, 17.0, 10.0);
    basicAtoms.back()->setId(0);
    
    basicAtoms.push_back(std::make_unique<Atom>("C"));
    basicAtoms.back()->setCoordinates(19.0, 16.2, 9.5);
    basicAtoms.back()->setId(1);
    
    basicAtoms.push_back(std::make_unique<Atom>("C"));
    basicAtoms.back()->setCoordinates(19.5, 15.0, 8.7);
    basicAtoms.back()->setId(2);
    
    basicAtoms.push_back(std::make_unique<Atom>("O"));
    basicAtoms.back()->setCoordinates(18.3, 14.5, 8.2);
    basicAtoms.back()->setId(3);
    
    basicAtoms.push_back(std::make_unique<Atom>("C"));
    basicAtoms.back()->setCoordinates(17.5, 13.0, 7.5);
    basicAtoms.back()->setId(4);
    
    // Enrich atoms with physical properties
    AtomBuilder builder;
    auto enrichedAtoms = builder.buildAtoms(basicAtoms);
    
    std::cout << "Created " << enrichedAtoms.size() << " atoms\n\n";
    
    // Test different voxel sizes
    std::vector<double> voxelSizes = {2.0, 1.0, 0.5};
    
    for (double voxelSize : voxelSizes) {
        std::cout << "\n--- Voxelization with voxel size = " << voxelSize << " Å ---\n";
        
        // Create voxel grid
        VoxelGrid voxelGrid(enrichedAtoms, voxelSize, 1.0);
        
        // Print statistics
        voxelGrid.printStatistics();
        
        // Generate hexahedral mesh
        HexMesh mesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);
        
        std::cout << "Generated hexahedral mesh:\n";
        std::cout << "  Nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "  Elements: " << mesh.getElementCount() << "\n";
        
        // Calculate node sharing efficiency
        int theoreticalNodes = voxelGrid.getOccupiedVoxelCount() * 8;
        double efficiency = 0.0;
        if (theoreticalNodes > 0) {
            efficiency = (1.0 - (double)mesh.getNodeCount() / theoreticalNodes) * 100.0;
        }
        std::cout << "  Node sharing efficiency: " << std::fixed << std::setprecision(1) 
                  << efficiency << "%\n";
        
        // Show first few occupied voxels
        const auto& occupiedVoxels = voxelGrid.getOccupiedVoxels();
        int numToShow = std::min(5, (int)occupiedVoxels.size());
        if (numToShow > 0) {
            std::cout << "\nFirst " << numToShow << " occupied voxels:\n";
            for (int i = 0; i < numToShow; ++i) {
                const Voxel& v = occupiedVoxels[i];
                std::cout << "  Voxel [" << v.index.i << "," << v.index.j << "," << v.index.k << "]: "
                          << "center=(" << std::fixed << std::setprecision(2)
                          << v.center.x << ", " << v.center.y << ", " << v.center.z << "), "
                          << "intersects " << v.atomIds.size() << " atom(s)\n";
            }
        }
    }
}

void demonstrateWithPDB() {
    std::cout << "\n\n=== Voxelization with PDB File ===\n\n";
    
    // Try to load test PDB file
    std::string pdbFile = "../data/test_peptide.pdb";
    
    try {
        // Parse PDB
        auto basicAtoms = PDBParser::parsePDBFile(pdbFile);
        std::cout << "Loaded " << basicAtoms.size() << " atoms from " << pdbFile << "\n\n";
        
        // Enrich with properties
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(basicAtoms);
        
        // Create voxel grid with 1.0 Å voxels
        double voxelSize = 1.0;
        VoxelGrid voxelGrid(enrichedAtoms, voxelSize, 2.0);
        
        voxelGrid.printStatistics();
        
        // Generate mesh
        HexMesh mesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);
        
        std::cout << "Generated hexahedral mesh:\n";
        std::cout << "  Nodes: " << mesh.getNodeCount() << "\n";
        std::cout << "  Elements: " << mesh.getElementCount() << "\n";
        
        // Calculate mesh statistics
        int theoreticalNodes = voxelGrid.getOccupiedVoxelCount() * 8;
        double efficiency = (1.0 - (double)mesh.getNodeCount() / theoreticalNodes) * 100.0;
        std::cout << "  Node sharing efficiency: " << std::fixed << std::setprecision(1) 
                  << efficiency << "%\n\n";
        
    } catch (const std::exception& e) {
        std::cout << "Could not load PDB file: " << e.what() << "\n";
        std::cout << "Skipping PDB demonstration.\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "BioMesh - Voxelization Demonstration\n";
    std::cout << "======================================\n\n";
    
    if (argc > 2) {
        // User provided PDB file and voxel size
        std::string pdbFile = argv[1];
        double voxelSize = std::stod(argv[2]);
        
        try {
            auto basicAtoms = PDBParser::parsePDBFile(pdbFile);
            std::cout << "Loaded " << basicAtoms.size() << " atoms from " << pdbFile << "\n";
            
            AtomBuilder builder;
            auto enrichedAtoms = builder.buildAtoms(basicAtoms);
            
            double padding = (argc > 3) ? std::stod(argv[3]) : 2.0;
            
            VoxelGrid voxelGrid(enrichedAtoms, voxelSize, padding);
            voxelGrid.printStatistics();
            
            HexMesh mesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);
            
            std::cout << "Hexahedral mesh:\n";
            std::cout << "  Nodes: " << mesh.getNodeCount() << "\n";
            std::cout << "  Elements: " << mesh.getElementCount() << "\n\n";
            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 1;
        }
        
    } else {
        // Run demonstrations
        demonstrateVoxelization();
        demonstrateWithPDB();
    }
    
    std::cout << "\nUsage: " << (argc > 0 ? argv[0] : "voxel_demo") 
              << " [pdb_file] [voxel_size] [padding]\n";
    std::cout << "  pdb_file: Path to PDB file\n";
    std::cout << "  voxel_size: Edge length of voxels in Angstroms\n";
    std::cout << "  padding: Optional padding around bounding box (default: 2.0 Å)\n\n";
    
    return 0;
}
