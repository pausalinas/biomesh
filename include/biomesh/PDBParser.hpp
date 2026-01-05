#pragma once

#include "biomesh/Atom.hpp"
#include <vector>
#include <string>
#include <memory>

namespace biomesh {

/**
 * @brief PDB file parser for extracting atom information
 * 
 * Reads PDB files and extracts atom coordinates and element information.
 */
class PDBParser {
public:
    /**
     * @brief Parse PDB file and extract atoms
     * @param filename Path to PDB file
     * @return Vector of parsed atoms (with coordinates but no physical properties)
     * @throws std::runtime_error if file cannot be opened or parsed
     */
    static std::vector<std::unique_ptr<Atom>> parsePDBFile(const std::string& filename);

    /**
     * @brief Parse PDB content from string
     * @param content PDB file content as string
     * @return Vector of parsed atoms
     * @throws std::runtime_error if content cannot be parsed
     */
    static std::vector<std::unique_ptr<Atom>> parsePDBContent(const std::string& content);

private:
    /**
     * @brief Parse a single ATOM line from PDB
     * @param line PDB ATOM line
     * @param atomId Sequential atom identifier
     * @return Parsed atom or nullptr if line is invalid
     */
    static std::unique_ptr<Atom> parseAtomLine(const std::string& line, size_t atomId);

    /**
     * @brief Extract element symbol from atom name with residue context
     * @param atomName Atom name from PDB file (columns 13-16)
     * @param residueName Residue name from PDB file (columns 18-20), optional
     * @return Chemical element symbol
     */
    static std::string extractElement(const std::string& atomName, const std::string& residueName = "");

    /**
     * @brief Validate element symbol against atomic database
     * @param element Element symbol to validate
     * @return true if element exists in atomic database
     */
    static bool isValidElement(const std::string& element);

    /**
     * @brief Resolve ambiguous atom names using residue context
     * @param atomName Atom name from PDB file
     * @param residueName Residue name from PDB file
     * @return Resolved chemical element symbol
     */
    static std::string resolveAmbiguousElement(const std::string& atomName, const std::string& residueName);

    /**
     * @brief Parse coordinate value from PDB line
     * @param line PDB line
     * @param start Start position
     * @param length Field length
     * @return Parsed coordinate value
     */
    static double parseCoordinate(const std::string& line, size_t start, size_t length);
};

} // namespace biomesh