#include "biomesh/AtomicSpec.hpp"
#include <stdexcept>

namespace biomesh {

AtomicSpecDatabase& AtomicSpecDatabase::getInstance() {
    static AtomicSpecDatabase instance;
    return instance;
}

AtomicSpecDatabase::AtomicSpecDatabase() {
    initializeDefaultSpecs();
}

void AtomicSpecDatabase::initializeDefaultSpecs() {
    // Van der Waals radii and atomic masses for common elements
    // Van der Waals radii data sources: Bondi (1964) and Rowland & Taylor (1996)
    // Atomic masses from NIST atomic masses (2020)
    
    specs_["H"] = AtomicSpec("H", 1.20, 1.008);    // Hydrogen
    specs_["C"] = AtomicSpec("C", 1.70, 12.011);   // Carbon
    specs_["N"] = AtomicSpec("N", 1.55, 14.007);   // Nitrogen
    specs_["O"] = AtomicSpec("O", 1.52, 15.999);   // Oxygen
    specs_["P"] = AtomicSpec("P", 1.80, 30.974);   // Phosphorus
    specs_["S"] = AtomicSpec("S", 1.80, 32.065);   // Sulfur
    specs_["F"] = AtomicSpec("F", 1.47, 18.998);   // Fluorine
    specs_["Cl"] = AtomicSpec("Cl", 1.75, 35.453); // Chlorine
    specs_["Br"] = AtomicSpec("Br", 1.85, 79.904); // Bromine
    specs_["I"] = AtomicSpec("I", 1.98, 126.904);  // Iodine
    specs_["Na"] = AtomicSpec("Na", 2.27, 22.990); // Sodium
    specs_["Mg"] = AtomicSpec("Mg", 1.73, 24.305); // Magnesium
    specs_["K"] = AtomicSpec("K", 2.75, 39.098);   // Potassium
    specs_["Ca"] = AtomicSpec("Ca", 2.31, 40.078); // Calcium
    specs_["Fe"] = AtomicSpec("Fe", 1.80, 55.845); // Iron
    specs_["Zn"] = AtomicSpec("Zn", 1.39, 65.38);  // Zinc
    specs_["Se"] = AtomicSpec("Se", 1.90, 78.96);  // Selenium
}

const AtomicSpec& AtomicSpecDatabase::getSpec(const std::string& element) const {
    auto it = specs_.find(element);
    if (it == specs_.end()) {
        throw std::runtime_error("Element '" + element + "' not found in atomic specification database");
    }
    return it->second;
}

bool AtomicSpecDatabase::hasElement(const std::string& element) const {
    return specs_.find(element) != specs_.end();
}

void AtomicSpecDatabase::addSpec(const AtomicSpec& spec) {
    specs_[spec.elementSymbol] = spec;
}

} // namespace biomesh