#pragma once

#include "biomesh/Atom.hpp"
#include "biomesh/ResidueClassifier.hpp"
#include <vector>
#include <memory>

namespace biomesh {

/**
 * @brief Filter type presets for common filtering scenarios
 */
enum class FilterType {
    All,              // Keep everything (default)
    ProteinOnly,      // Only amino acids
    NucleicAcidOnly,  // Only DNA/RNA
    NoWater,          // Remove water molecules
    NoHeteroatoms,    // Remove all HETATM records
    Custom            // User-defined filter
};

/**
 * @brief Filter for selecting atoms based on molecule type
 * 
 * Provides preset filters and a builder pattern for custom filtering.
 * Filters atoms based on residue classification (proteins, nucleic acids,
 * water, ions, etc.).
 */
class MoleculeFilter {
public:
    /**
     * @brief Default constructor - keeps everything
     */
    MoleculeFilter();

    // Preset static factory methods
    
    /**
     * @brief Create filter that keeps all atoms (no filtering)
     * @return Filter that keeps everything
     */
    static MoleculeFilter all();
    
    /**
     * @brief Create filter that keeps only protein atoms
     * @return Filter that keeps only amino acids
     */
    static MoleculeFilter proteinOnly();
    
    /**
     * @brief Create filter that keeps only nucleic acid atoms
     * @return Filter that keeps only DNA/RNA
     */
    static MoleculeFilter nucleicAcidOnly();
    
    /**
     * @brief Create filter that removes water molecules
     * @return Filter that keeps everything except water
     */
    static MoleculeFilter noWater();

    // Builder pattern for custom filters
    
    /**
     * @brief Set whether to keep protein atoms
     * @param keep true to keep proteins, false to filter them out
     * @return Reference to this filter for chaining
     */
    MoleculeFilter& setKeepProteins(bool keep);
    
    /**
     * @brief Set whether to keep nucleic acid atoms
     * @param keep true to keep nucleic acids, false to filter them out
     * @return Reference to this filter for chaining
     */
    MoleculeFilter& setKeepNucleicAcids(bool keep);
    
    /**
     * @brief Set whether to keep water molecules
     * @param keep true to keep water, false to filter it out
     * @return Reference to this filter for chaining
     */
    MoleculeFilter& setKeepWater(bool keep);
    
    /**
     * @brief Set whether to keep ion atoms
     * @param keep true to keep ions, false to filter them out
     * @return Reference to this filter for chaining
     */
    MoleculeFilter& setKeepIons(bool keep);
    
    /**
     * @brief Set whether to keep other/unknown residues
     * @param keep true to keep others, false to filter them out
     * @return Reference to this filter for chaining
     */
    MoleculeFilter& setKeepOthers(bool keep);

    // Main filtering methods
    
    /**
     * @brief Filter a collection of atoms
     * @param atoms Input atoms to filter
     * @return Vector of atoms that pass the filter (deep copies)
     */
    std::vector<std::unique_ptr<Atom>> filter(const std::vector<std::unique_ptr<Atom>>& atoms) const;
    
    /**
     * @brief Check if a single atom should be kept
     * @param atom Atom to check
     * @return true if atom passes the filter
     */
    bool shouldKeep(const Atom& atom) const;

private:
    FilterType type_;
    bool keepProteins_;
    bool keepNucleicAcids_;
    bool keepWater_;
    bool keepIons_;
    bool keepOthers_;
};

} // namespace biomesh
