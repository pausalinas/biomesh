#pragma once

#include <string>
#include <set>

namespace biomesh {

/**
 * @brief Classifier for identifying molecule types based on residue names
 * 
 * Provides static methods to classify residues as proteins, nucleic acids,
 * water, ions, and other biomolecule types commonly found in PDB files.
 */
class ResidueClassifier {
public:
    /**
     * @brief Check if residue is a standard or non-standard amino acid (protein)
     * @param residueName Three-letter residue name from PDB file
     * @return true if residue is an amino acid
     */
    static bool isProtein(const std::string& residueName);

    /**
     * @brief Check if residue is a nucleic acid (DNA or RNA)
     * @param residueName Residue name from PDB file
     * @return true if residue is DNA or RNA
     */
    static bool isNucleicAcid(const std::string& residueName);

    /**
     * @brief Check if residue is DNA
     * @param residueName Residue name from PDB file
     * @return true if residue is DNA
     */
    static bool isDNA(const std::string& residueName);

    /**
     * @brief Check if residue is RNA
     * @param residueName Residue name from PDB file
     * @return true if residue is RNA
     */
    static bool isRNA(const std::string& residueName);

    /**
     * @brief Check if residue is water
     * @param residueName Residue name from PDB file
     * @return true if residue is water
     */
    static bool isWater(const std::string& residueName);

    /**
     * @brief Check if residue is a common ion
     * @param residueName Residue name from PDB file
     * @return true if residue is an ion
     */
    static bool isIon(const std::string& residueName);

private:
    // Standard amino acids (20 common + non-standard)
    static const std::set<std::string> aminoAcids_;
    
    // DNA residues
    static const std::set<std::string> dnaResidues_;
    
    // RNA residues
    static const std::set<std::string> rnaResidues_;
    
    // Water molecules
    static const std::set<std::string> waterResidues_;
    
    // Common ions
    static const std::set<std::string> ionResidues_;
};

} // namespace biomesh
