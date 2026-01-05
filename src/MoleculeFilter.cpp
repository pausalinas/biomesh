#include "biomesh/MoleculeFilter.hpp"

namespace biomesh {

MoleculeFilter::MoleculeFilter() 
    : type_(FilterType::All),
      keepProteins_(true),
      keepNucleicAcids_(true),
      keepWater_(true),
      keepIons_(true),
      keepOthers_(true) {
}

MoleculeFilter MoleculeFilter::all() {
    MoleculeFilter filter;
    filter.type_ = FilterType::All;
    filter.keepProteins_ = true;
    filter.keepNucleicAcids_ = true;
    filter.keepWater_ = true;
    filter.keepIons_ = true;
    filter.keepOthers_ = true;
    return filter;
}

MoleculeFilter MoleculeFilter::proteinOnly() {
    MoleculeFilter filter;
    filter.type_ = FilterType::ProteinOnly;
    filter.keepProteins_ = true;
    filter.keepNucleicAcids_ = false;
    filter.keepWater_ = false;
    filter.keepIons_ = false;
    filter.keepOthers_ = false;
    return filter;
}

MoleculeFilter MoleculeFilter::nucleicAcidOnly() {
    MoleculeFilter filter;
    filter.type_ = FilterType::NucleicAcidOnly;
    filter.keepProteins_ = false;
    filter.keepNucleicAcids_ = true;
    filter.keepWater_ = false;
    filter.keepIons_ = false;
    filter.keepOthers_ = false;
    return filter;
}

MoleculeFilter MoleculeFilter::noWater() {
    MoleculeFilter filter;
    filter.type_ = FilterType::NoWater;
    filter.keepProteins_ = true;
    filter.keepNucleicAcids_ = true;
    filter.keepWater_ = false;
    filter.keepIons_ = true;
    filter.keepOthers_ = true;
    return filter;
}

MoleculeFilter& MoleculeFilter::setKeepProteins(bool keep) {
    keepProteins_ = keep;
    type_ = FilterType::Custom;
    return *this;
}

MoleculeFilter& MoleculeFilter::setKeepNucleicAcids(bool keep) {
    keepNucleicAcids_ = keep;
    type_ = FilterType::Custom;
    return *this;
}

MoleculeFilter& MoleculeFilter::setKeepWater(bool keep) {
    keepWater_ = keep;
    type_ = FilterType::Custom;
    return *this;
}

MoleculeFilter& MoleculeFilter::setKeepIons(bool keep) {
    keepIons_ = keep;
    type_ = FilterType::Custom;
    return *this;
}

MoleculeFilter& MoleculeFilter::setKeepOthers(bool keep) {
    keepOthers_ = keep;
    type_ = FilterType::Custom;
    return *this;
}

std::vector<std::unique_ptr<Atom>> MoleculeFilter::filter(const std::vector<std::unique_ptr<Atom>>& atoms) const {
    std::vector<std::unique_ptr<Atom>> filtered;
    
    for (const auto& atom : atoms) {
        if (shouldKeep(*atom)) {
            // Create a deep copy of the atom
            auto newAtom = std::make_unique<Atom>(
                atom->getChemicalElement(),
                atom->getAtomicRadius(),
                atom->getAtomicMass()
            );
            newAtom->setCoordinates(atom->getX(), atom->getY(), atom->getZ());
            newAtom->setId(atom->getId());
            newAtom->setResidueName(atom->getResidueName());
            newAtom->setAtomName(atom->getAtomName());
            newAtom->setResidueNumber(atom->getResidueNumber());
            newAtom->setChainID(atom->getChainID());
            
            filtered.push_back(std::move(newAtom));
        }
    }
    
    return filtered;
}

bool MoleculeFilter::shouldKeep(const Atom& atom) const {
    const std::string& residueName = atom.getResidueName();
    
    // Check each molecule type
    if (ResidueClassifier::isProtein(residueName)) {
        return keepProteins_;
    }
    
    if (ResidueClassifier::isNucleicAcid(residueName)) {
        return keepNucleicAcids_;
    }
    
    if (ResidueClassifier::isWater(residueName)) {
        return keepWater_;
    }
    
    if (ResidueClassifier::isIon(residueName)) {
        return keepIons_;
    }
    
    // If not classified as any known type, treat as "other"
    return keepOthers_;
}

} // namespace biomesh
