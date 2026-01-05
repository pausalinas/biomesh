#include "biomesh/BioMesh.hpp"
#include <iostream>
#include <iomanip>

using namespace biomesh;

/**
 * Example demonstrating the complete workflow with MoleculeFilter
 * 
 * Shows how to integrate filtering into the BioMesh pipeline:
 * PDB File → Parse → Filter → Enrich → Calculate Bounding Box
 */
int main(int argc, char* argv[]) {
    std::cout << "=== BioMesh Workflow with MoleculeFilter ===" << std::endl << std::endl;
    
    // Use command line argument or default test file
    std::string pdbFile = "../data/mixed_molecules.pdb";
    if (argc > 1) {
        pdbFile = argv[1];
    }
    
    try {
        // Step 1: Parse PDB file (extracts all atoms)
        std::cout << "Step 1: Parsing PDB file..." << std::endl;
        auto atoms = PDBParser::parsePDBFile(pdbFile);
        std::cout << "  Parsed " << atoms.size() << " atoms" << std::endl << std::endl;
        
        // Step 2: Filter molecules (optional preprocessing step)
        std::cout << "Step 2: Filtering molecules..." << std::endl;
        std::cout << "  Keeping: proteins and nucleic acids" << std::endl;
        std::cout << "  Removing: water, ions, and other molecules" << std::endl;
        
        MoleculeFilter filter;
        filter.setKeepProteins(true)
              .setKeepNucleicAcids(true)
              .setKeepWater(false)
              .setKeepIons(false)
              .setKeepOthers(false);
        
        auto filteredAtoms = filter.filter(atoms);
        std::cout << "  Filtered to " << filteredAtoms.size() << " atoms" << std::endl << std::endl;
        
        // Step 3: Enrich atoms with physical properties
        std::cout << "Step 3: Enriching atoms with physical properties..." << std::endl;
        AtomBuilder builder;
        auto enrichedAtoms = builder.buildAtoms(filteredAtoms);
        std::cout << "  Enriched " << enrichedAtoms.size() << " atoms" << std::endl << std::endl;
        
        // Step 4: Calculate bounding box
        std::cout << "Step 4: Calculating bounding box..." << std::endl;
        double padding = 2.0;  // 2 Angstrom padding
        BoundingBox bbox(enrichedAtoms, padding);
        std::cout << "  Bounding box calculated with " << padding << " Å padding" << std::endl << std::endl;
        
        // Display results
        std::cout << "=== Results ===" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        
        Point3D min = bbox.getMin();
        Point3D max = bbox.getMax();
        Point3D center = bbox.getCenter();
        Point3D dims = bbox.getDimensions();
        
        std::cout << "Bounding Box:" << std::endl;
        std::cout << "  Min corner: (" << min.x << ", " << min.y << ", " << min.z << ")" << std::endl;
        std::cout << "  Max corner: (" << max.x << ", " << max.y << ", " << max.z << ")" << std::endl;
        std::cout << "  Center:     (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
        std::cout << "  Dimensions: " << dims.x << " × " << dims.y << " × " << dims.z << " Å" << std::endl;
        std::cout << "  Volume:     " << bbox.getVolume() << " Ų" << std::endl;
        std::cout << "  Surface:    " << bbox.getSurfaceArea() << " Ų" << std::endl << std::endl;
        
        // Show atom count breakdown by type
        std::cout << "Atom Count by Type:" << std::endl;
        int proteinCount = 0, nucleicCount = 0;
        for (const auto& atom : enrichedAtoms) {
            if (ResidueClassifier::isProtein(atom->getResidueName())) {
                proteinCount++;
            } else if (ResidueClassifier::isNucleicAcid(atom->getResidueName())) {
                nucleicCount++;
            }
        }
        std::cout << "  Protein atoms:      " << proteinCount << std::endl;
        std::cout << "  Nucleic acid atoms: " << nucleicCount << std::endl;
        std::cout << "  Total:              " << enrichedAtoms.size() << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n✓ Workflow completed successfully!" << std::endl;
    return 0;
}
