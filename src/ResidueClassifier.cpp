#include "biomesh/ResidueClassifier.hpp"
#include <algorithm>

namespace biomesh {

// Standard amino acids (20 common + non-standard MSE, SEC, PYL)
const std::set<std::string> ResidueClassifier::aminoAcids_ = {
    // Standard 20 amino acids
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
    "THR", "TRP", "TYR", "VAL",
    // Non-standard amino acids
    "MSE",  // Selenomethionine
    "SEC",  // Selenocysteine
    "PYL"   // Pyrrolysine
};

// DNA residues (deoxy)
const std::set<std::string> ResidueClassifier::dnaResidues_ = {
    "DA",   // Deoxyadenosine
    "DT",   // Deoxythymidine
    "DG",   // Deoxyguanosine
    "DC",   // Deoxycytidine
    "DU"    // Deoxyuridine
};

// RNA residues
const std::set<std::string> ResidueClassifier::rnaResidues_ = {
    "A",    // Adenosine
    "U",    // Uridine
    "G",    // Guanosine
    "C",    // Cytidine
    "ADE",  // Adenosine (alternative)
    "URA",  // Uracil (alternative)
    "GUA",  // Guanine (alternative)
    "CYT"   // Cytosine (alternative)
};

// Water molecules
const std::set<std::string> ResidueClassifier::waterResidues_ = {
    "HOH",   // Standard water in PDB
    "WAT",   // Water
    "H2O",   // Water
    "SOL",   // Solvent (often water)
    "TIP",   // TIP water model
    "TIP3",  // TIP3P water model
    "TIP4"   // TIP4P water model
};

// Common ions
const std::set<std::string> ResidueClassifier::ionResidues_ = {
    "NA",   // Sodium
    "CL",   // Chloride
    "K",    // Potassium
    "CA",   // Calcium
    "MG",   // Magnesium
    "ZN",   // Zinc
    "FE",   // Iron
    "CU",   // Copper
    "MN"    // Manganese
};

bool ResidueClassifier::isProtein(const std::string& residueName) {
    // Convert to uppercase for case-insensitive comparison
    std::string upperName = residueName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    return aminoAcids_.find(upperName) != aminoAcids_.end();
}

bool ResidueClassifier::isNucleicAcid(const std::string& residueName) {
    return isDNA(residueName) || isRNA(residueName);
}

bool ResidueClassifier::isDNA(const std::string& residueName) {
    std::string upperName = residueName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    return dnaResidues_.find(upperName) != dnaResidues_.end();
}

bool ResidueClassifier::isRNA(const std::string& residueName) {
    std::string upperName = residueName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    return rnaResidues_.find(upperName) != rnaResidues_.end();
}

bool ResidueClassifier::isWater(const std::string& residueName) {
    std::string upperName = residueName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    return waterResidues_.find(upperName) != waterResidues_.end();
}

bool ResidueClassifier::isIon(const std::string& residueName) {
    std::string upperName = residueName;
    std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
    
    return ionResidues_.find(upperName) != ionResidues_.end();
}

} // namespace biomesh
