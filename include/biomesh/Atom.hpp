#pragma once

#include <string>
#include <memory>

namespace biomesh {

/**
 * @brief Enhanced Atom structure with physical properties
 * 
 * Represents an atom with chemical element, coordinates, and physical properties.
 * Uses van der Waals radii for accurate molecular representation.
 */
class Atom {
public:
    /**
     * @brief Constructor with chemical element only
     * @param element Chemical element symbol (e.g., "C", "N", "O")
     */
    explicit Atom(const std::string& element);

    /**
     * @brief Constructor with chemical element and radius
     * @param element Chemical element symbol
     * @param radius Van der Waals radius in Angstroms
     */
    Atom(const std::string& element, double radius);

    /**
     * @brief Constructor with chemical element, radius, and atomic mass
     * @param element Chemical element symbol
     * @param radius Van der Waals radius in Angstroms
     * @param mass Atomic mass in Daltons
     */
    Atom(const std::string& element, double radius, double mass);

    // Getters
    const std::string& getChemicalElement() const { return chemicalElement_; }
    double getX() const { return x_; }
    double getY() const { return y_; }
    double getZ() const { return z_; }
    double getAtomicRadius() const { return atomicRadius_; }
    double getAtomicMass() const { return atomicMass_; }
    size_t getId() const { return id_; }
    const std::string& getResidueName() const { return residueName_; }
    const std::string& getAtomName() const { return atomName_; }
    int getResidueNumber() const { return residueNumber_; }
    char getChainID() const { return chainID_; }

    // Setters for coordinates
    void setCoordinates(double x, double y, double z);
    void setId(size_t id) { id_ = id; }
    void setResidueName(const std::string& residueName) { residueName_ = residueName; }
    void setAtomName(const std::string& atomName) { atomName_ = atomName; }
    void setResidueNumber(int residueNumber) { residueNumber_ = residueNumber; }
    void setChainID(char chainID) { chainID_ = chainID; }

private:
    std::string chemicalElement_;
    double x_ = 0.0;
    double y_ = 0.0;
    double z_ = 0.0;
    double atomicRadius_ = 0.0;
    double atomicMass_ = 0.0;
    size_t id_ = 0;
    std::string residueName_;
    std::string atomName_;
    int residueNumber_ = 0;
    char chainID_ = ' ';
};

} // namespace biomesh