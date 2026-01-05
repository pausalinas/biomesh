#include "biomesh/AtomBuilder.hpp"
#include <stdexcept>
#include <set>

namespace biomesh {

std::vector<std::unique_ptr<Atom>> AtomBuilder::buildAtoms(
    const std::vector<std::unique_ptr<Atom>>& basicAtoms) const {
    
    std::vector<std::unique_ptr<Atom>> enrichedAtoms;
    enrichedAtoms.reserve(basicAtoms.size());

    for (const auto& basicAtom : basicAtoms) {
        enrichedAtoms.push_back(buildAtom(*basicAtom));
    }

    return enrichedAtoms;
}

std::unique_ptr<Atom> AtomBuilder::buildAtom(const Atom& basicAtom) const {
    const std::string& element = basicAtom.getChemicalElement();
    
    if (!database_.hasElement(element)) {
        throw std::runtime_error("Element '" + element + "' not found in atomic specification database");
    }

    const AtomicSpec& spec = database_.getSpec(element);
    
    // Create new atom with full properties
    auto enrichedAtom = std::make_unique<Atom>(element, spec.radius, spec.mass);
    enrichedAtom->setCoordinates(basicAtom.getX(), basicAtom.getY(), basicAtom.getZ());
    enrichedAtom->setId(basicAtom.getId());
    
    // Preserve residue information
    enrichedAtom->setResidueName(basicAtom.getResidueName());
    enrichedAtom->setAtomName(basicAtom.getAtomName());
    enrichedAtom->setResidueNumber(basicAtom.getResidueNumber());
    enrichedAtom->setChainID(basicAtom.getChainID());

    return enrichedAtom;
}

bool AtomBuilder::areAllElementsSupported(const std::vector<std::unique_ptr<Atom>>& atoms) const {
    for (const auto& atom : atoms) {
        if (!database_.hasElement(atom->getChemicalElement())) {
            return false;
        }
    }
    return true;
}

std::vector<std::string> AtomBuilder::getUnsupportedElements(
    const std::vector<std::unique_ptr<Atom>>& atoms) const {
    
    std::set<std::string> unsupportedSet;
    
    for (const auto& atom : atoms) {
        const std::string& element = atom->getChemicalElement();
        if (!database_.hasElement(element)) {
            unsupportedSet.insert(element);
        }
    }

    return std::vector<std::string>(unsupportedSet.begin(), unsupportedSet.end());
}

} // namespace biomesh