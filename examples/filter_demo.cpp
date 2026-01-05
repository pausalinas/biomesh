#include "biomesh/PDBParser.hpp"
#include "biomesh/MoleculeFilter.hpp"
#include "biomesh/ResidueClassifier.hpp"
#include <iostream>
#include <iomanip>

using namespace biomesh;

int main(int argc, char* argv[]) {
    std::cout << "=== MoleculeFilter Demonstration ===" << std::endl << std::endl;
    
    // Parse PDB file
    std::string pdbFile = "../data/test_peptide.pdb";
    if (argc > 1) {
        pdbFile = argv[1];
    }
    
    std::cout << "Parsing PDB file: " << pdbFile << std::endl;
    
    try {
        auto atoms = PDBParser::parsePDBFile(pdbFile);
        
        std::cout << "Total atoms parsed: " << atoms.size() << std::endl << std::endl;
        
        // Display atom information
        std::cout << "=== Atom Details ===" << std::endl;
        std::cout << std::setw(5) << "ID" 
                  << std::setw(8) << "Element"
                  << std::setw(10) << "ResName"
                  << std::setw(10) << "AtomName"
                  << std::setw(8) << "ResNum"
                  << std::setw(8) << "Chain"
                  << std::endl;
        std::cout << std::string(55, '-') << std::endl;
        
        for (const auto& atom : atoms) {
            std::cout << std::setw(5) << atom->getId()
                      << std::setw(8) << atom->getChemicalElement()
                      << std::setw(10) << atom->getResidueName()
                      << std::setw(10) << atom->getAtomName()
                      << std::setw(8) << atom->getResidueNumber()
                      << std::setw(8) << atom->getChainID()
                      << std::endl;
        }
        std::cout << std::endl;
        
        // Test filter presets
        std::cout << "=== Filter Tests ===" << std::endl;
        
        // Test 1: All (no filtering)
        auto filterAll = MoleculeFilter::all();
        auto filteredAll = filterAll.filter(atoms);
        std::cout << "All atoms: " << filteredAll.size() << std::endl;
        
        // Test 2: Protein only
        auto filterProtein = MoleculeFilter::proteinOnly();
        auto filteredProtein = filterProtein.filter(atoms);
        std::cout << "Protein only: " << filteredProtein.size() << std::endl;
        
        // Test 3: No water
        auto filterNoWater = MoleculeFilter::noWater();
        auto filteredNoWater = filterNoWater.filter(atoms);
        std::cout << "No water: " << filteredNoWater.size() << std::endl;
        
        // Test 4: Custom filter
        MoleculeFilter customFilter;
        customFilter.setKeepProteins(true)
                    .setKeepWater(false)
                    .setKeepIons(false);
        auto filteredCustom = customFilter.filter(atoms);
        std::cout << "Custom (proteins, no water/ions): " << filteredCustom.size() << std::endl;
        std::cout << std::endl;
        
        // Classify residues
        std::cout << "=== Residue Classification ===" << std::endl;
        int proteinCount = 0, nucleicCount = 0, waterCount = 0, ionCount = 0, otherCount = 0;
        
        for (const auto& atom : atoms) {
            const std::string& res = atom->getResidueName();
            if (ResidueClassifier::isProtein(res)) {
                proteinCount++;
            } else if (ResidueClassifier::isNucleicAcid(res)) {
                nucleicCount++;
            } else if (ResidueClassifier::isWater(res)) {
                waterCount++;
            } else if (ResidueClassifier::isIon(res)) {
                ionCount++;
            } else {
                otherCount++;
            }
        }
        
        std::cout << "Protein atoms: " << proteinCount << std::endl;
        std::cout << "Nucleic acid atoms: " << nucleicCount << std::endl;
        std::cout << "Water atoms: " << waterCount << std::endl;
        std::cout << "Ion atoms: " << ionCount << std::endl;
        std::cout << "Other atoms: " << otherCount << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
