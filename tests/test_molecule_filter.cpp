#include <gtest/gtest.h>
#include "biomesh/ResidueClassifier.hpp"
#include "biomesh/MoleculeFilter.hpp"
#include "biomesh/Atom.hpp"
#include <memory>
#include <vector>

using namespace biomesh;

// Test fixtures
class ResidueClassifierTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class MoleculeFilterTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test atoms with different residue types
        
        // Protein atoms
        auto protein1 = std::make_unique<Atom>("C", 1.7, 12.011);
        protein1->setCoordinates(1.0, 2.0, 3.0);
        protein1->setResidueName("ALA");
        protein1->setAtomName("CA");
        testAtoms.push_back(std::move(protein1));
        
        auto protein2 = std::make_unique<Atom>("N", 1.55, 14.007);
        protein2->setCoordinates(2.0, 3.0, 4.0);
        protein2->setResidueName("GLY");
        protein2->setAtomName("N");
        testAtoms.push_back(std::move(protein2));
        
        // DNA atoms
        auto dna1 = std::make_unique<Atom>("C", 1.7, 12.011);
        dna1->setCoordinates(3.0, 4.0, 5.0);
        dna1->setResidueName("DA");
        dna1->setAtomName("C1'");
        testAtoms.push_back(std::move(dna1));
        
        // RNA atoms
        auto rna1 = std::make_unique<Atom>("C", 1.7, 12.011);
        rna1->setCoordinates(4.0, 5.0, 6.0);
        rna1->setResidueName("A");
        rna1->setAtomName("C1'");
        testAtoms.push_back(std::move(rna1));
        
        // Water molecules
        auto water1 = std::make_unique<Atom>("O", 1.52, 15.999);
        water1->setCoordinates(5.0, 6.0, 7.0);
        water1->setResidueName("HOH");
        water1->setAtomName("O");
        testAtoms.push_back(std::move(water1));
        
        auto water2 = std::make_unique<Atom>("O", 1.52, 15.999);
        water2->setCoordinates(6.0, 7.0, 8.0);
        water2->setResidueName("WAT");
        water2->setAtomName("O");
        testAtoms.push_back(std::move(water2));
        
        // Ions
        auto ion1 = std::make_unique<Atom>("Na", 2.27, 22.990);
        ion1->setCoordinates(7.0, 8.0, 9.0);
        ion1->setResidueName("NA");
        ion1->setAtomName("NA");
        testAtoms.push_back(std::move(ion1));
        
        auto ion2 = std::make_unique<Atom>("Ca", 2.31, 40.078);
        ion2->setCoordinates(8.0, 9.0, 10.0);
        ion2->setResidueName("CA");
        ion2->setAtomName("CA");
        testAtoms.push_back(std::move(ion2));
        
        // Unknown/other residue
        auto other1 = std::make_unique<Atom>("C", 1.7, 12.011);
        other1->setCoordinates(9.0, 10.0, 11.0);
        other1->setResidueName("UNK");
        other1->setAtomName("C1");
        testAtoms.push_back(std::move(other1));
    }
    
    void TearDown() override {
        testAtoms.clear();
    }
    
    std::vector<std::unique_ptr<Atom>> testAtoms;
};

// ResidueClassifier Tests

TEST_F(ResidueClassifierTest, IdentifiesStandardAminoAcids) {
    EXPECT_TRUE(ResidueClassifier::isProtein("ALA"));
    EXPECT_TRUE(ResidueClassifier::isProtein("ARG"));
    EXPECT_TRUE(ResidueClassifier::isProtein("ASN"));
    EXPECT_TRUE(ResidueClassifier::isProtein("ASP"));
    EXPECT_TRUE(ResidueClassifier::isProtein("CYS"));
    EXPECT_TRUE(ResidueClassifier::isProtein("GLN"));
    EXPECT_TRUE(ResidueClassifier::isProtein("GLU"));
    EXPECT_TRUE(ResidueClassifier::isProtein("GLY"));
    EXPECT_TRUE(ResidueClassifier::isProtein("HIS"));
    EXPECT_TRUE(ResidueClassifier::isProtein("ILE"));
    EXPECT_TRUE(ResidueClassifier::isProtein("LEU"));
    EXPECT_TRUE(ResidueClassifier::isProtein("LYS"));
    EXPECT_TRUE(ResidueClassifier::isProtein("MET"));
    EXPECT_TRUE(ResidueClassifier::isProtein("PHE"));
    EXPECT_TRUE(ResidueClassifier::isProtein("PRO"));
    EXPECT_TRUE(ResidueClassifier::isProtein("SER"));
    EXPECT_TRUE(ResidueClassifier::isProtein("THR"));
    EXPECT_TRUE(ResidueClassifier::isProtein("TRP"));
    EXPECT_TRUE(ResidueClassifier::isProtein("TYR"));
    EXPECT_TRUE(ResidueClassifier::isProtein("VAL"));
}

TEST_F(ResidueClassifierTest, IdentifiesNonStandardAminoAcids) {
    EXPECT_TRUE(ResidueClassifier::isProtein("MSE"));  // Selenomethionine
    EXPECT_TRUE(ResidueClassifier::isProtein("SEC"));  // Selenocysteine
    EXPECT_TRUE(ResidueClassifier::isProtein("PYL"));  // Pyrrolysine
}

TEST_F(ResidueClassifierTest, CaseInsensitiveAminoAcids) {
    EXPECT_TRUE(ResidueClassifier::isProtein("ala"));
    EXPECT_TRUE(ResidueClassifier::isProtein("Ala"));
    EXPECT_TRUE(ResidueClassifier::isProtein("ALA"));
}

TEST_F(ResidueClassifierTest, IdentifiesDNA) {
    EXPECT_TRUE(ResidueClassifier::isDNA("DA"));
    EXPECT_TRUE(ResidueClassifier::isDNA("DT"));
    EXPECT_TRUE(ResidueClassifier::isDNA("DG"));
    EXPECT_TRUE(ResidueClassifier::isDNA("DC"));
    EXPECT_TRUE(ResidueClassifier::isDNA("DU"));
}

TEST_F(ResidueClassifierTest, IdentifiesRNA) {
    EXPECT_TRUE(ResidueClassifier::isRNA("A"));
    EXPECT_TRUE(ResidueClassifier::isRNA("U"));
    EXPECT_TRUE(ResidueClassifier::isRNA("G"));
    EXPECT_TRUE(ResidueClassifier::isRNA("C"));
    EXPECT_TRUE(ResidueClassifier::isRNA("ADE"));
    EXPECT_TRUE(ResidueClassifier::isRNA("URA"));
    EXPECT_TRUE(ResidueClassifier::isRNA("GUA"));
    EXPECT_TRUE(ResidueClassifier::isRNA("CYT"));
}

TEST_F(ResidueClassifierTest, IdentifiesNucleicAcids) {
    EXPECT_TRUE(ResidueClassifier::isNucleicAcid("DA"));
    EXPECT_TRUE(ResidueClassifier::isNucleicAcid("A"));
    EXPECT_FALSE(ResidueClassifier::isNucleicAcid("ALA"));
}

TEST_F(ResidueClassifierTest, IdentifiesWater) {
    EXPECT_TRUE(ResidueClassifier::isWater("HOH"));
    EXPECT_TRUE(ResidueClassifier::isWater("WAT"));
    EXPECT_TRUE(ResidueClassifier::isWater("H2O"));
    EXPECT_TRUE(ResidueClassifier::isWater("SOL"));
    EXPECT_TRUE(ResidueClassifier::isWater("TIP"));
    EXPECT_TRUE(ResidueClassifier::isWater("TIP3"));
    EXPECT_TRUE(ResidueClassifier::isWater("TIP4"));
}

TEST_F(ResidueClassifierTest, CaseInsensitiveWater) {
    EXPECT_TRUE(ResidueClassifier::isWater("hoh"));
    EXPECT_TRUE(ResidueClassifier::isWater("Hoh"));
}

TEST_F(ResidueClassifierTest, IdentifiesIons) {
    EXPECT_TRUE(ResidueClassifier::isIon("NA"));
    EXPECT_TRUE(ResidueClassifier::isIon("CL"));
    EXPECT_TRUE(ResidueClassifier::isIon("K"));
    EXPECT_TRUE(ResidueClassifier::isIon("CA"));
    EXPECT_TRUE(ResidueClassifier::isIon("MG"));
    EXPECT_TRUE(ResidueClassifier::isIon("ZN"));
    EXPECT_TRUE(ResidueClassifier::isIon("FE"));
    EXPECT_TRUE(ResidueClassifier::isIon("CU"));
    EXPECT_TRUE(ResidueClassifier::isIon("MN"));
}

TEST_F(ResidueClassifierTest, RejectsUnknownResidues) {
    EXPECT_FALSE(ResidueClassifier::isProtein("XXX"));
    EXPECT_FALSE(ResidueClassifier::isNucleicAcid("XXX"));
    EXPECT_FALSE(ResidueClassifier::isWater("XXX"));
    EXPECT_FALSE(ResidueClassifier::isIon("XXX"));
}

// MoleculeFilter Tests

TEST_F(MoleculeFilterTest, FilterAllKeepsEverything) {
    auto filter = MoleculeFilter::all();
    auto filtered = filter.filter(testAtoms);
    
    EXPECT_EQ(testAtoms.size(), filtered.size());
    EXPECT_EQ(9, filtered.size());
}

TEST_F(MoleculeFilterTest, FilterProteinOnly) {
    auto filter = MoleculeFilter::proteinOnly();
    auto filtered = filter.filter(testAtoms);
    
    // Should keep 2 protein atoms (ALA, GLY)
    EXPECT_EQ(2, filtered.size());
    
    for (const auto& atom : filtered) {
        EXPECT_TRUE(ResidueClassifier::isProtein(atom->getResidueName()));
    }
}

TEST_F(MoleculeFilterTest, FilterNucleicAcidOnly) {
    auto filter = MoleculeFilter::nucleicAcidOnly();
    auto filtered = filter.filter(testAtoms);
    
    // Should keep 2 nucleic acid atoms (DA, A)
    EXPECT_EQ(2, filtered.size());
    
    for (const auto& atom : filtered) {
        EXPECT_TRUE(ResidueClassifier::isNucleicAcid(atom->getResidueName()));
    }
}

TEST_F(MoleculeFilterTest, FilterNoWater) {
    auto filter = MoleculeFilter::noWater();
    auto filtered = filter.filter(testAtoms);
    
    // Should remove 2 water molecules, keep 7 others
    EXPECT_EQ(7, filtered.size());
    
    for (const auto& atom : filtered) {
        EXPECT_FALSE(ResidueClassifier::isWater(atom->getResidueName()));
    }
}

TEST_F(MoleculeFilterTest, CustomFilterProteinAndDNA) {
    MoleculeFilter filter;
    filter.setKeepProteins(true)
          .setKeepNucleicAcids(true)
          .setKeepWater(false)
          .setKeepIons(false)
          .setKeepOthers(false);
    
    auto filtered = filter.filter(testAtoms);
    
    // Should keep 2 protein + 2 nucleic acid = 4 atoms
    EXPECT_EQ(4, filtered.size());
}

TEST_F(MoleculeFilterTest, CustomFilterNoIons) {
    MoleculeFilter filter;
    filter.setKeepProteins(true)
          .setKeepNucleicAcids(true)
          .setKeepWater(true)
          .setKeepIons(false)
          .setKeepOthers(true);
    
    auto filtered = filter.filter(testAtoms);
    
    // Should remove 2 ions, keep 7 others
    EXPECT_EQ(7, filtered.size());
    
    for (const auto& atom : filtered) {
        EXPECT_FALSE(ResidueClassifier::isIon(atom->getResidueName()));
    }
}

TEST_F(MoleculeFilterTest, ShouldKeepMethod) {
    auto filter = MoleculeFilter::proteinOnly();
    
    // Test protein atom
    Atom proteinAtom("C");
    proteinAtom.setResidueName("ALA");
    EXPECT_TRUE(filter.shouldKeep(proteinAtom));
    
    // Test water atom
    Atom waterAtom("O");
    waterAtom.setResidueName("HOH");
    EXPECT_FALSE(filter.shouldKeep(waterAtom));
}

TEST_F(MoleculeFilterTest, FilterPreservesAtomProperties) {
    auto filter = MoleculeFilter::proteinOnly();
    auto filtered = filter.filter(testAtoms);
    
    ASSERT_GT(filtered.size(), 0);
    
    // Check that properties are preserved
    auto& atom = filtered[0];
    EXPECT_EQ("C", atom->getChemicalElement());
    EXPECT_EQ(1.7, atom->getAtomicRadius());
    EXPECT_EQ(12.011, atom->getAtomicMass());
    EXPECT_EQ(1.0, atom->getX());
    EXPECT_EQ(2.0, atom->getY());
    EXPECT_EQ(3.0, atom->getZ());
    EXPECT_EQ("ALA", atom->getResidueName());
    EXPECT_EQ("CA", atom->getAtomName());
}

TEST_F(MoleculeFilterTest, EmptyInputReturnsEmpty) {
    std::vector<std::unique_ptr<Atom>> emptyAtoms;
    auto filter = MoleculeFilter::proteinOnly();
    auto filtered = filter.filter(emptyAtoms);
    
    EXPECT_EQ(0, filtered.size());
}

TEST_F(MoleculeFilterTest, UnknownResiduesHandledByOthers) {
    MoleculeFilter filter;
    filter.setKeepProteins(false)
          .setKeepNucleicAcids(false)
          .setKeepWater(false)
          .setKeepIons(false)
          .setKeepOthers(true);
    
    auto filtered = filter.filter(testAtoms);
    
    // Should keep only the 1 unknown residue (UNK)
    EXPECT_EQ(1, filtered.size());
    EXPECT_EQ("UNK", filtered[0]->getResidueName());
}

TEST_F(MoleculeFilterTest, BuilderPatternChaining) {
    // Test that builder pattern works with method chaining
    MoleculeFilter filter;
    auto& result = filter.setKeepProteins(true)
                         .setKeepWater(false)
                         .setKeepIons(false);
    
    // Should return reference to same object
    EXPECT_EQ(&filter, &result);
}
