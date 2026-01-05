#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomicSpec.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>
#include <vector>

namespace biomesh {

std::vector<std::unique_ptr<Atom>> PDBParser::parsePDBFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open PDB file: " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    file.close();

    return parsePDBContent(buffer.str());
}

std::vector<std::unique_ptr<Atom>> PDBParser::parsePDBContent(const std::string& content) {
    std::vector<std::unique_ptr<Atom>> atoms;
    std::istringstream stream(content);
    std::string line;
    size_t atomId = 0;

    while (std::getline(stream, line)) {
        // Check if line starts with "ATOM" or "HETATM"
        if (line.length() >= 4 && (line.substr(0, 4) == "ATOM" || 
            (line.length() >= 6 && line.substr(0, 6) == "HETATM"))) {
            auto atom = parseAtomLine(line, atomId++);
            if (atom) {
                atoms.push_back(std::move(atom));
            }
        }
    }

    if (atoms.empty()) {
        throw std::runtime_error("No valid ATOM or HETATM records found in PDB content");
    }

    return atoms;
}

std::unique_ptr<Atom> PDBParser::parseAtomLine(const std::string& line, size_t atomId) {
    // PDB format specification for ATOM records:
    // Columns  Data type      Field        Definition
    // 1-6      Record name    "ATOM  "
    // 7-11     Integer        serial       Atom serial number
    // 13-16    Atom           name         Atom name
    // 17       Character      altLoc       Alternate location indicator
    // 18-20    Residue name   resName      Residue name
    // 22       Character      chainID      Chain identifier
    // 23-26    Integer        resSeq       Residue sequence number
    // 31-38    Real(8.3)      x            Orthogonal coordinates for X in Angstroms
    // 39-46    Real(8.3)      y            Orthogonal coordinates for Y in Angstroms
    // 47-54    Real(8.3)      z            Orthogonal coordinates for Z in Angstroms
    // 77-78    LString(2)     element      Element symbol

    if (line.length() < 54) {
        return nullptr; // Line too short to contain coordinates
    }

    try {
        // Extract coordinates
        double x = parseCoordinate(line, 30, 8);  // positions 31-38 (0-indexed: 30-37)
        double y = parseCoordinate(line, 38, 8);  // positions 39-46 (0-indexed: 38-45)
        double z = parseCoordinate(line, 46, 8);  // positions 47-54 (0-indexed: 46-53)

        // Extract atom name (columns 13-16)
        std::string atomName;
        if (line.length() >= 16) {
            atomName = line.substr(12, 4);
            // Remove leading/trailing whitespace but preserve internal spaces
            atomName.erase(0, atomName.find_first_not_of(' '));
            atomName.erase(atomName.find_last_not_of(' ') + 1);
        }

        // Extract residue name for context-aware element detection
        std::string residueName;
        if (line.length() >= 20) {
            residueName = line.substr(17, 3);  // positions 18-20 (0-indexed: 17-19)
            // Remove whitespace
            residueName.erase(std::remove_if(residueName.begin(), residueName.end(), ::isspace), residueName.end());
        }

        // Extract chain ID (column 22)
        char chainID = ' ';
        if (line.length() >= 22) {
            chainID = line[21];  // position 22 (0-indexed: 21)
        }

        // Extract residue number (columns 23-26)
        int residueNumber = 0;
        if (line.length() >= 26) {
            std::string resNumStr = line.substr(22, 4);  // positions 23-26 (0-indexed: 22-25)
            resNumStr.erase(std::remove_if(resNumStr.begin(), resNumStr.end(), ::isspace), resNumStr.end());
            if (!resNumStr.empty()) {
                try {
                    residueNumber = std::stoi(resNumStr);
                } catch (...) {
                    // Keep as 0 if parsing fails
                }
            }
        }

        // Extract element symbol - try columns 77-78 first (PDB standard)
        std::string element;
        if (line.length() >= 78) {
            // Try to get element from columns 77-78
            element = line.substr(76, 2);
            // Remove whitespace
            element.erase(std::remove_if(element.begin(), element.end(), ::isspace), element.end());
            
            // Validate the element against atomic database
            if (!element.empty() && !isValidElement(element)) {
                element.clear(); // Invalid element, fall back to atom name parsing
            }
        }
        
        // Fallback: extract from atom name with residue context
        if (element.empty() && !atomName.empty()) {
            element = extractElement(atomName, residueName);
        }

        if (element.empty()) {
            return nullptr; // Could not determine element
        }

        // Create atom with just the element (no properties yet)
        auto atom = std::make_unique<Atom>(element);
        atom->setCoordinates(x, y, z);
        atom->setId(atomId);
        atom->setResidueName(residueName);
        atom->setAtomName(atomName);
        atom->setResidueNumber(residueNumber);
        atom->setChainID(chainID);

        return atom;
    } catch (const std::exception&) {
        return nullptr; // Parsing failed
    }
}

std::string PDBParser::extractElement(const std::string& atomName, const std::string& residueName) {
    // Remove whitespace
    std::string name = atomName;
    name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());
    
    if (name.empty()) {
        return "";
    }

    // First, try to resolve ambiguous cases using residue context
    std::string resolvedElement = resolveAmbiguousElement(name, residueName);
    if (!resolvedElement.empty() && isValidElement(resolvedElement)) {
        return resolvedElement;
    }

    // Standard element extraction logic
    std::string element;
    element += std::toupper(name[0]);
    
    // Check for two-letter elements
    if (name.length() > 1) {
        char second = std::tolower(name[1]);
        
        // Common two-letter elements in biological systems
        std::string twoLetterElement = element + second;
        if (isValidElement(twoLetterElement)) {
            return twoLetterElement;
        }
    }

    // Validate single-letter element
    if (isValidElement(element)) {
        return element;
    }

    return ""; // Could not determine valid element
}

bool PDBParser::isValidElement(const std::string& element) {
    return AtomicSpecDatabase::getInstance().hasElement(element);
}

std::string PDBParser::resolveAmbiguousElement(const std::string& atomName, const std::string& residueName) {
    // Handle common ambiguous cases with residue context
    
    // CA can be alpha carbon (in amino acids) or calcium (in metal sites)
    if (atomName == "CA") {
        // Check if we're in a standard amino acid residue
        static const std::vector<std::string> aminoAcids = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", 
            "THR", "TRP", "TYR", "VAL", "SEC", "PYL"
        };
        
        for (const auto& aa : aminoAcids) {
            if (residueName == aa) {
                return "C"; // Alpha carbon in amino acid
            }
        }
        
        // If not in amino acid, likely calcium
        return "Ca";
    }
    
    // Handle other potential ambiguities
    // (Add more cases as needed based on common PDB patterns)
    
    return ""; // No specific resolution
}

double PDBParser::parseCoordinate(const std::string& line, size_t start, size_t length) {
    if (start + length > line.length()) {
        throw std::runtime_error("Line too short for coordinate extraction");
    }

    std::string coordStr = line.substr(start, length);
    
    // Remove leading/trailing whitespace
    coordStr.erase(0, coordStr.find_first_not_of(" \t"));
    coordStr.erase(coordStr.find_last_not_of(" \t") + 1);
    
    if (coordStr.empty()) {
        throw std::runtime_error("Empty coordinate field");
    }

    return std::stod(coordStr);
}

} // namespace biomesh