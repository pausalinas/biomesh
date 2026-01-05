#pragma once

/**
 * @file BioMesh.hpp
 * @brief Main header file for the BioMesh C++ module
 * 
 * This module provides functionality for:
 * - Parsing PDB structure files
 * - Extracting and enriching atom information with physical properties
 * - Filtering molecules by type (proteins, nucleic acids, water, ions)
 * - Computing molecular bounding boxes
 * 
 * @author BioMesh Team
 * @version 1.0.0
 */

#include "biomesh/Atom.hpp"
#include "biomesh/AtomicSpec.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomBuilder.hpp"
#include "biomesh/BoundingBox.hpp"
#include "biomesh/ResidueClassifier.hpp"
#include "biomesh/MoleculeFilter.hpp"

namespace biomesh {

/**
 * @brief Main workflow function for processing PDB files
 * 
 * This convenience function combines all the steps:
 * 1. Parse PDB file
 * 2. Enrich atoms with physical properties
 * 3. Calculate bounding box
 * 
 * @param pdbFilename Path to PDB file
 * @param padding Additional padding for bounding box (default: 0.0)
 * @return Pair of (enriched atoms, bounding box)
 * @throws std::runtime_error if file cannot be processed
 */
std::pair<std::vector<std::unique_ptr<Atom>>, BoundingBox> 
processPDBFile(const std::string& pdbFilename, double padding = 0.0);

} // namespace biomesh