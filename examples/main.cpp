#include "biomesh/BioMesh.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using namespace biomesh;

void printAtomInfo(const std::vector<std::unique_ptr<Atom>>& atoms) {
    std::cout << "\n=== Atom Information ===" << std::endl;
    std::cout << std::setw(6) << "ID" 
              << std::setw(8) << "Element"
              << std::setw(12) << "X (Å)"
              << std::setw(12) << "Y (Å)"
              << std::setw(12) << "Z (Å)"
              << std::setw(12) << "Radius (Å)"
              << std::setw(12) << "Mass (Da)" << std::endl;
    std::cout << std::string(74, '-') << std::endl;

    for (const auto& atom : atoms) {
        std::cout << std::setw(6) << atom->getId()
                  << std::setw(8) << atom->getChemicalElement()
                  << std::setw(12) << std::fixed << std::setprecision(3) << atom->getX()
                  << std::setw(12) << std::fixed << std::setprecision(3) << atom->getY()
                  << std::setw(12) << std::fixed << std::setprecision(3) << atom->getZ()
                  << std::setw(12) << std::fixed << std::setprecision(3) << atom->getAtomicRadius()
                  << std::setw(12) << std::fixed << std::setprecision(3) << atom->getAtomicMass()
                  << std::endl;
    }
}

void printBoundingBoxInfo(const BoundingBox& bbox) {
    std::cout << "\n=== Bounding Box Information ===" << std::endl;
    
    Point3D min = bbox.getMin();
    Point3D max = bbox.getMax();
    Point3D center = bbox.getCenter();
    Point3D dimensions = bbox.getDimensions();
    
    std::cout << "Min corner: (" << std::fixed << std::setprecision(3) 
              << min.x << ", " << min.y << ", " << min.z << ")" << std::endl;
    std::cout << "Max corner: (" << std::fixed << std::setprecision(3) 
              << max.x << ", " << max.y << ", " << max.z << ")" << std::endl;
    std::cout << "Center:     (" << std::fixed << std::setprecision(3) 
              << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
    std::cout << "Dimensions: " << std::fixed << std::setprecision(3) 
              << dimensions.x << " × " << dimensions.y << " × " << dimensions.z << " Å" << std::endl;
    std::cout << "Volume:     " << std::fixed << std::setprecision(3) 
              << bbox.getVolume() << " Ų" << std::endl;
    std::cout << "Surface:    " << std::fixed << std::setprecision(3) 
              << bbox.getSurfaceArea() << " Ų" << std::endl;
}

void demonstrateStepByStep() {
    std::cout << "\n=== Step-by-Step Demonstration ===" << std::endl;
    
    // Create a simple test PDB content
    std::string testPdbContent = 
        "HEADER    TEST MOLECULE                           01-JAN-70   TEST            \n"
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
        "ATOM      3  C   ALA A   1      18.500  15.300  10.600  1.00 20.00           C  \n"
        "ATOM      4  O   ALA A   1      17.400  14.800  10.500  1.00 20.00           O  \n"
        "ATOM      5  CB  ALA A   1      17.900  17.100   8.900  1.00 20.00           C  \n"
        "END                                                                             \n";

    try {
        // Step 1: Parse PDB content
        std::cout << "\nStep 1: Parsing PDB content..." << std::endl;
        auto basicAtoms = PDBParser::parsePDBContent(testPdbContent);
        std::cout << "Parsed " << basicAtoms.size() << " atoms" << std::endl;

        // Step 2: Enrich atoms with properties
        std::cout << "\nStep 2: Enriching atoms with physical properties..." << std::endl;
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(basicAtoms);
        std::cout << "Enriched " << enrichedAtoms.size() << " atoms" << std::endl;

        // Step 3: Calculate bounding box
        std::cout << "\nStep 3: Calculating bounding box..." << std::endl;
        double padding = 2.0; // 2 Å padding
        BoundingBox boundingBox(enrichedAtoms, padding);
        std::cout << "Bounding box calculated with " << padding << " Å padding" << std::endl;

        // Display results
        printAtomInfo(enrichedAtoms);
        printBoundingBoxInfo(boundingBox);

    } catch (const std::exception& e) {
        std::cerr << "Error in step-by-step demonstration: " << e.what() << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << "BioMesh C++ Module - PDB Parser and Bounding Box Calculator" << std::endl;
    std::cout << "===========================================================" << std::endl;

    // Run step-by-step demonstration with built-in test data
    demonstrateStepByStep();

    // If a PDB file is provided as argument, process it
    if (argc > 1) {
        std::string pdbFile = argv[1];
        double padding = (argc > 2) ? std::stod(argv[2]) : 0.0;

        std::cout << "\n\n=== Processing PDB File: " << pdbFile << " ===" << std::endl;
        
        try {
            auto [atoms, boundingBox] = processPDBFile(pdbFile, padding);
            
            std::cout << "Successfully processed " << atoms.size() << " atoms" << std::endl;
            
            // For large files, show only first few atoms
            if (atoms.size() > 10) {
                std::vector<std::unique_ptr<Atom>> firstTen;
                for (size_t i = 0; i < 10; ++i) {
                    // Create a copy for display (atoms vector is moved)
                    const auto& atom = atoms[i];
                    auto displayAtom = std::make_unique<Atom>(
                        atom->getChemicalElement(), 
                        atom->getAtomicRadius(), 
                        atom->getAtomicMass()
                    );
                    displayAtom->setCoordinates(atom->getX(), atom->getY(), atom->getZ());
                    displayAtom->setId(atom->getId());
                    firstTen.push_back(std::move(displayAtom));
                }
                std::cout << "\nShowing first 10 atoms:" << std::endl;
                printAtomInfo(firstTen);
                std::cout << "\n... (" << (atoms.size() - 10) << " more atoms)" << std::endl;
            } else {
                printAtomInfo(atoms);
            }
            
            printBoundingBoxInfo(boundingBox);
            
        } catch (const std::exception& e) {
            std::cerr << "Error processing PDB file: " << e.what() << std::endl;
            return 1;
        }
    } else {
        std::cout << "\nUsage: " << argv[0] << " [pdb_file] [padding]" << std::endl;
        std::cout << "  pdb_file: Path to PDB file (optional)" << std::endl;
        std::cout << "  padding:  Additional padding for bounding box in Angstroms (optional, default: 0.0)" << std::endl;
    }

    return 0;
}