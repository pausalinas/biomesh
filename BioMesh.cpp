#include "biomesh/BioMesh.hpp"

namespace biomesh {

std::pair<std::vector<std::unique_ptr<Atom>>, BoundingBox> 
processPDBFile(const std::string& pdbFilename, double padding) {
    // Step 1: Parse PDB file
    auto basicAtoms = PDBParser::parsePDBFile(pdbFilename);
    
    // Step 2: Enrich atoms with physical properties
    AtomBuilder builder;
    auto enrichedAtoms = builder.buildAtoms(basicAtoms);
    
    // Step 3: Calculate bounding box
    BoundingBox boundingBox(enrichedAtoms, padding);
    
    return std::make_pair(std::move(enrichedAtoms), boundingBox);
}

} // namespace biomesh