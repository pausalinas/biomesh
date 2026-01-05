#pragma once

#include "biomesh/Atom.hpp"
#include "biomesh/AtomicSpec.hpp"
#include <vector>
#include <memory>

namespace biomesh {

/**
 * @brief Builder class for enriching atoms with physical properties
 * 
 * Takes atoms with basic information and enriches them with van der Waals
 * radius and mass from the atomic specification database.
 * Uses van der Waals radii for accurate molecular representation.
 */
class AtomBuilder {
public:
    /**
     * @brief Default constructor
     */
    AtomBuilder() = default;

    /**
     * @brief Build fully initialized atoms from basic atoms
     * @param basicAtoms Vector of atoms with element and coordinates only
     * @return Vector of fully initialized atoms with physical properties
     * @throws std::runtime_error if any element is not found in the database
     */
    std::vector<std::unique_ptr<Atom>> buildAtoms(
        const std::vector<std::unique_ptr<Atom>>& basicAtoms) const;

    /**
     * @brief Build a single fully initialized atom
     * @param basicAtom Atom with element and coordinates only
     * @return Fully initialized atom with physical properties
     * @throws std::runtime_error if element is not found in the database
     */
    std::unique_ptr<Atom> buildAtom(const Atom& basicAtom) const;

    /**
     * @brief Check if all elements in the atom vector are supported
     * @param atoms Vector of atoms to check
     * @return true if all elements are in the database
     */
    bool areAllElementsSupported(const std::vector<std::unique_ptr<Atom>>& atoms) const;

    /**
     * @brief Get list of unsupported elements
     * @param atoms Vector of atoms to check
     * @return Vector of unsupported element symbols
     */
    std::vector<std::string> getUnsupportedElements(
        const std::vector<std::unique_ptr<Atom>>& atoms) const;

private:
    AtomicSpecDatabase& database_ = AtomicSpecDatabase::getInstance();
};

} // namespace biomesh