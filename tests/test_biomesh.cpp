#include <gtest/gtest.h>
#include "biomesh/Atom.hpp"
#include "biomesh/AtomicSpec.hpp"
#include "biomesh/AtomBuilder.hpp"
#include "biomesh/BoundingBox.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/VoxelMeshGenerator.hpp"
#include <memory>
#include <vector>
#include <set>

using namespace biomesh;

// Test fixtures
class AtomTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class AtomicSpecTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class AtomBuilderTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class BoundingBoxTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class PDBParserTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class VoxelGridTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

class VoxelMeshGeneratorTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

// Atom class tests
TEST_F(AtomTest, ConstructorWithElementOnly) {
    Atom atom("C");
    EXPECT_EQ("C", atom.getChemicalElement());
    EXPECT_EQ(0.0, atom.getX());
    EXPECT_EQ(0.0, atom.getY());
    EXPECT_EQ(0.0, atom.getZ());
    EXPECT_EQ(0.0, atom.getAtomicRadius());
    EXPECT_EQ(0.0, atom.getAtomicMass());
    EXPECT_EQ(0, atom.getId());
}

TEST_F(AtomTest, ConstructorWithElementAndRadius) {
    Atom atom("N", 1.55);
    EXPECT_EQ("N", atom.getChemicalElement());
    EXPECT_EQ(1.55, atom.getAtomicRadius());
    EXPECT_EQ(0.0, atom.getAtomicMass());
}

TEST_F(AtomTest, ConstructorWithAllProperties) {
    Atom atom("O", 1.52, 15.999);
    EXPECT_EQ("O", atom.getChemicalElement());
    EXPECT_EQ(1.52, atom.getAtomicRadius());
    EXPECT_EQ(15.999, atom.getAtomicMass());
}

TEST_F(AtomTest, SetCoordinates) {
    Atom atom("H");
    atom.setCoordinates(1.0, 2.0, 3.0);
    EXPECT_EQ(1.0, atom.getX());
    EXPECT_EQ(2.0, atom.getY());
    EXPECT_EQ(3.0, atom.getZ());
}

TEST_F(AtomTest, SetId) {
    Atom atom("C");
    atom.setId(42);
    EXPECT_EQ(42, atom.getId());
}

// AtomicSpec tests
TEST_F(AtomicSpecTest, DatabaseSingleton) {
    auto& db1 = AtomicSpecDatabase::getInstance();
    auto& db2 = AtomicSpecDatabase::getInstance();
    EXPECT_EQ(&db1, &db2);
}

TEST_F(AtomicSpecTest, DefaultSpecsExist) {
    auto& db = AtomicSpecDatabase::getInstance();
    
    // Test common elements
    EXPECT_TRUE(db.hasElement("H"));
    EXPECT_TRUE(db.hasElement("C"));
    EXPECT_TRUE(db.hasElement("N"));
    EXPECT_TRUE(db.hasElement("O"));
    EXPECT_TRUE(db.hasElement("P"));
    EXPECT_TRUE(db.hasElement("S"));
    
    // Test some specific values - using van der Waals radii
    const auto& carbon = db.getSpec("C");
    EXPECT_EQ("C", carbon.elementSymbol);
    EXPECT_EQ(1.70, carbon.radius);  // Van der Waals radius
    EXPECT_EQ(12.011, carbon.mass);
    
    const auto& hydrogen = db.getSpec("H");
    EXPECT_EQ("H", hydrogen.elementSymbol);
    EXPECT_EQ(1.20, hydrogen.radius);  // Van der Waals radius
    EXPECT_EQ(1.008, hydrogen.mass);
}

TEST_F(AtomicSpecTest, UnknownElementThrows) {
    auto& db = AtomicSpecDatabase::getInstance();
    EXPECT_FALSE(db.hasElement("Xx"));
    EXPECT_THROW(db.getSpec("Xx"), std::runtime_error);
}

TEST_F(AtomicSpecTest, AddCustomSpec) {
    auto& db = AtomicSpecDatabase::getInstance();
    AtomicSpec custom("X", 2.0, 100.0);
    db.addSpec(custom);
    
    EXPECT_TRUE(db.hasElement("X"));
    const auto& retrieved = db.getSpec("X");
    EXPECT_EQ("X", retrieved.elementSymbol);
    EXPECT_EQ(2.0, retrieved.radius);
    EXPECT_EQ(100.0, retrieved.mass);
}

// AtomBuilder tests - focusing on property assignment
TEST_F(AtomBuilderTest, CorrectPropertyAssignment) {
    // Create basic atoms
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    auto carbon = std::make_unique<Atom>("C");
    carbon->setCoordinates(1.0, 2.0, 3.0);
    carbon->setId(1);
    
    auto nitrogen = std::make_unique<Atom>("N");
    nitrogen->setCoordinates(4.0, 5.0, 6.0);
    nitrogen->setId(2);
    
    basicAtoms.push_back(std::move(carbon));
    basicAtoms.push_back(std::move(nitrogen));
    
    // Build enriched atoms
    AtomBuilder builder;
    auto enrichedAtoms = builder.buildAtoms(basicAtoms);
    
    ASSERT_EQ(2, enrichedAtoms.size());
    
    // Check carbon properties - using van der Waals radii
    const auto& enrichedCarbon = enrichedAtoms[0];
    EXPECT_EQ("C", enrichedCarbon->getChemicalElement());
    EXPECT_EQ(1.0, enrichedCarbon->getX());
    EXPECT_EQ(2.0, enrichedCarbon->getY());
    EXPECT_EQ(3.0, enrichedCarbon->getZ());
    EXPECT_EQ(1.70, enrichedCarbon->getAtomicRadius());  // Van der Waals radius
    EXPECT_EQ(12.011, enrichedCarbon->getAtomicMass());
    EXPECT_EQ(1, enrichedCarbon->getId());
    
    // Check nitrogen properties - using van der Waals radii
    const auto& enrichedNitrogen = enrichedAtoms[1];
    EXPECT_EQ("N", enrichedNitrogen->getChemicalElement());
    EXPECT_EQ(4.0, enrichedNitrogen->getX());
    EXPECT_EQ(5.0, enrichedNitrogen->getY());
    EXPECT_EQ(6.0, enrichedNitrogen->getZ());
    EXPECT_EQ(1.55, enrichedNitrogen->getAtomicRadius());  // Van der Waals radius
    EXPECT_EQ(14.007, enrichedNitrogen->getAtomicMass());
    EXPECT_EQ(2, enrichedNitrogen->getId());
}

TEST_F(AtomBuilderTest, UnsupportedElementThrows) {
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    auto unknownAtom = std::make_unique<Atom>("Xx");
    basicAtoms.push_back(std::move(unknownAtom));
    
    AtomBuilder builder;
    EXPECT_THROW(builder.buildAtoms(basicAtoms), std::runtime_error);
}

TEST_F(AtomBuilderTest, UnsupportedElementDetection) {
    std::vector<std::unique_ptr<Atom>> basicAtoms;
    basicAtoms.push_back(std::make_unique<Atom>("C"));
    basicAtoms.push_back(std::make_unique<Atom>("Xx"));
    basicAtoms.push_back(std::make_unique<Atom>("Yy"));
    
    AtomBuilder builder;
    EXPECT_FALSE(builder.areAllElementsSupported(basicAtoms));
    
    auto unsupported = builder.getUnsupportedElements(basicAtoms);
    EXPECT_EQ(2, unsupported.size());
    EXPECT_TRUE(std::find(unsupported.begin(), unsupported.end(), "Xx") != unsupported.end());
    EXPECT_TRUE(std::find(unsupported.begin(), unsupported.end(), "Yy") != unsupported.end());
}

// BoundingBox tests - focusing on correct calculation
TEST_F(BoundingBoxTest, CorrectBoundingBoxCalculation) {
    // Create test atoms with known coordinates and radii
    std::vector<std::unique_ptr<Atom>> atoms;
    
    // Atom at origin with radius 1.0
    auto atom1 = std::make_unique<Atom>("C", 1.0, 12.0);
    atom1->setCoordinates(0.0, 0.0, 0.0);
    atoms.push_back(std::move(atom1));
    
    // Atom at (10, 0, 0) with radius 2.0
    auto atom2 = std::make_unique<Atom>("N", 2.0, 14.0);
    atom2->setCoordinates(10.0, 0.0, 0.0);
    atoms.push_back(std::move(atom2));
    
    // No padding
    BoundingBox bbox(atoms, 0.0);
    
    // Expected bounds:
    // Atom1 at (0,0,0) with radius 1.0: bounds (-1,-1,-1) to (1,1,1)
    // Atom2 at (10,0,0) with radius 2.0: bounds (8,-2,-2) to (12,2,2)
    // Combined: Min: (-1, -2, -2), Max: (12, 2, 2)
    Point3D expectedMin(-1.0, -2.0, -2.0);
    Point3D expectedMax(12.0, 2.0, 2.0);
    
    EXPECT_NEAR(expectedMin.x, bbox.getMin().x, 1e-6);
    EXPECT_NEAR(expectedMin.y, bbox.getMin().y, 1e-6);
    EXPECT_NEAR(expectedMin.z, bbox.getMin().z, 1e-6);
    EXPECT_NEAR(expectedMax.x, bbox.getMax().x, 1e-6);
    EXPECT_NEAR(expectedMax.y, bbox.getMax().y, 1e-6);
    EXPECT_NEAR(expectedMax.z, bbox.getMax().z, 1e-6);
}

TEST_F(BoundingBoxTest, BoundingBoxWithPadding) {
    std::vector<std::unique_ptr<Atom>> atoms;
    auto atom = std::make_unique<Atom>("C", 1.0, 12.0);
    atom->setCoordinates(0.0, 0.0, 0.0);
    atoms.push_back(std::move(atom));
    
    double padding = 2.0;
    BoundingBox bbox(atoms, padding);
    
    // Expected bounds with padding:
    // Min: (-1-2, -1-2, -1-2) = (-3, -3, -3)
    // Max: (1+2, 1+2, 1+2) = (3, 3, 3)
    EXPECT_NEAR(-3.0, bbox.getMin().x, 1e-6);
    EXPECT_NEAR(-3.0, bbox.getMin().y, 1e-6);
    EXPECT_NEAR(-3.0, bbox.getMin().z, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().x, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().y, 1e-6);
    EXPECT_NEAR(3.0, bbox.getMax().z, 1e-6);
}

TEST_F(BoundingBoxTest, BoundingBoxProperties) {
    std::vector<std::unique_ptr<Atom>> atoms;
    auto atom = std::make_unique<Atom>("C", 0.0, 12.0); // No radius for simple calculation
    atom->setCoordinates(1.0, 2.0, 3.0);
    atoms.push_back(std::move(atom));
    
    BoundingBox bbox(atoms, 1.0); // 1.0 padding
    
    // Expected bounds: (0, 1, 2) to (2, 3, 4)
    Point3D expectedCenter(1.0, 2.0, 3.0);
    Point3D expectedDimensions(2.0, 2.0, 2.0);
    double expectedVolume = 8.0;
    double expectedSurfaceArea = 24.0;
    
    Point3D center = bbox.getCenter();
    Point3D dimensions = bbox.getDimensions();
    
    EXPECT_NEAR(expectedCenter.x, center.x, 1e-6);
    EXPECT_NEAR(expectedCenter.y, center.y, 1e-6);
    EXPECT_NEAR(expectedCenter.z, center.z, 1e-6);
    EXPECT_NEAR(expectedDimensions.x, dimensions.x, 1e-6);
    EXPECT_NEAR(expectedDimensions.y, dimensions.y, 1e-6);
    EXPECT_NEAR(expectedDimensions.z, dimensions.z, 1e-6);
    EXPECT_NEAR(expectedVolume, bbox.getVolume(), 1e-6);
    EXPECT_NEAR(expectedSurfaceArea, bbox.getSurfaceArea(), 1e-6);
}

TEST_F(BoundingBoxTest, EmptyAtomVectorThrows) {
    std::vector<std::unique_ptr<Atom>> emptyAtoms;
    EXPECT_THROW(BoundingBox bbox(emptyAtoms), std::runtime_error);
}

// PDB Parser tests
TEST_F(PDBParserTest, ParseSimplePDBContent) {
    std::string pdbContent = 
        "HEADER    TEST                                    01-JAN-70   TEST            \n"
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
        "END                                                                             \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(2, atoms.size());
    
    // Check first atom
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_NEAR(20.154, atoms[0]->getX(), 1e-3);
    EXPECT_NEAR(16.967, atoms[0]->getY(), 1e-3);
    EXPECT_NEAR(10.000, atoms[0]->getZ(), 1e-3);
    EXPECT_EQ(0, atoms[0]->getId());
    
    // Check second atom
    EXPECT_EQ("C", atoms[1]->getChemicalElement());
    EXPECT_NEAR(19.030, atoms[1]->getX(), 1e-3);
    EXPECT_NEAR(16.200, atoms[1]->getY(), 1e-3);
    EXPECT_NEAR(9.500, atoms[1]->getZ(), 1e-3);
    EXPECT_EQ(1, atoms[1]->getId());
}

TEST_F(PDBParserTest, EmptyContentThrows) {
    std::string emptyContent = "HEADER    TEST\nEND\n";
    EXPECT_THROW(PDBParser::parsePDBContent(emptyContent), std::runtime_error);
}

TEST_F(PDBParserTest, NonexistentFileThrows) {
    EXPECT_THROW(PDBParser::parsePDBFile("nonexistent_file.pdb"), std::runtime_error);
}

// Enhanced PDB Parser tests for improved element extraction
TEST_F(PDBParserTest, ElementExtractionFromColumns77_78) {
    // Test element extraction when columns 77-78 are present and valid
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00           C  \n"
        "ATOM      3  FE  HEM A   2      18.000  15.000   8.000  1.00 30.00          Fe  \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_EQ("C", atoms[1]->getChemicalElement());
    EXPECT_EQ("Fe", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, ElementExtractionMissingColumns77_78) {
    // Test element extraction when columns 77-78 are missing (shorter lines)
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00\n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00\n"
        "ATOM      3  O   ALA A   1      17.500  14.500   7.000  1.00 20.00\n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());
    EXPECT_EQ("C", atoms[1]->getChemicalElement());  // CA -> C (alpha carbon in amino acid)
    EXPECT_EQ("O", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, AmbiguousCAResolution) {
    // Test CA disambiguation: CA in amino acid vs Ca (calcium)
    std::string pdbContentAminoAcid = 
        "ATOM      1  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00\n";
    
    std::string pdbContentCalcium = 
        "ATOM      1  CA  CAL A   1      19.030  16.200   9.500  1.00 20.00\n";
    
    auto atomsAA = PDBParser::parsePDBContent(pdbContentAminoAcid);
    auto atomsCa = PDBParser::parsePDBContent(pdbContentCalcium);
    
    ASSERT_EQ(1, atomsAA.size());
    ASSERT_EQ(1, atomsCa.size());
    EXPECT_EQ("C", atomsAA[0]->getChemicalElement());   // Alpha carbon in alanine
    EXPECT_EQ("Ca", atomsCa[0]->getChemicalElement());  // Calcium in non-amino acid
}

TEST_F(PDBParserTest, TwoLetterElementExtraction) {
    // Test extraction of two-letter elements
    std::string pdbContent = 
        "ATOM      1  MG  HEM A   1      20.000  16.000  10.000  1.00 20.00\n"
        "ATOM      2  ZN  ZNC A   2      19.000  15.000   9.000  1.00 20.00\n"
        "ATOM      3  CL  CLA A   3      18.000  14.000   8.000  1.00 20.00\n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(3, atoms.size());
    EXPECT_EQ("Mg", atoms[0]->getChemicalElement());
    EXPECT_EQ("Zn", atoms[1]->getChemicalElement());
    EXPECT_EQ("Cl", atoms[2]->getChemicalElement());
}

TEST_F(PDBParserTest, InvalidElementFallback) {
    // Test fallback when columns 77-78 contain invalid elements
    std::string pdbContent = 
        "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00          XX  \n"
        "ATOM      2  CA  ALA A   1      19.030  16.200   9.500  1.00 20.00          YY  \n";
    
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    
    ASSERT_EQ(2, atoms.size());
    EXPECT_EQ("N", atoms[0]->getChemicalElement());   // Fallback to atom name
    EXPECT_EQ("C", atoms[1]->getChemicalElement());   // Fallback to atom name (CA -> C in amino acid)
}

TEST_F(PDBParserTest, ElementValidationAgainstDatabase) {
    // Test that only elements in the atomic database are accepted
    // Use 'Q' which is not a real element and shouldn't be in the database
    std::string pdbContent = 
        "ATOM      1  QQ  UNK A   1      20.000  16.000  10.000  1.00 20.00\n";
    
    // This should throw because no valid atoms can be parsed (QQ -> Q, and Q is not in our database)
    EXPECT_THROW(PDBParser::parsePDBContent(pdbContent), std::runtime_error);
}

TEST_F(PDBParserTest, MixedValidInvalidElements) {
    // Test with mixed valid and invalid elements
    // Use 'Q' which is not a real element and shouldn't be in the database
    std::string pdbContent = 
        "ATOM      1  C   ALA A   1      20.000  16.000  10.000  1.00 20.00\n"
        "ATOM      2  QQ  UNK A   2      21.000  17.000  11.000  1.00 20.00\n";
    
    // Should parse only the valid atom (carbon)
    auto atoms = PDBParser::parsePDBContent(pdbContent);
    ASSERT_EQ(1, atoms.size());
    EXPECT_EQ("C", atoms[0]->getChemicalElement());
}

// VoxelGrid tests
TEST_F(VoxelGridTest, BasicGridCreation) {
    // Create test atoms
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    // Create voxel grid
    VoxelGrid grid(atoms, 1.0, 0.5);
    
    // Check grid properties
    EXPECT_GT(grid.getTotalVoxelCount(), 0);
    EXPECT_GT(grid.getOccupiedVoxelCount(), 0);
    EXPECT_EQ(grid.getVoxelSize(), 1.0);
}

TEST_F(VoxelGridTest, VoxelOccupancy) {
    // Create a carbon atom at origin with radius ~0.67 Å
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(5.0, 5.0, 5.0);
    atoms.back()->setId(0);
    
    // Create grid with 1 Å voxels
    VoxelGrid grid(atoms, 1.0, 1.0);
    
    // Check that at least one voxel is occupied
    EXPECT_GT(grid.getOccupiedVoxelCount(), 0);
    
    // Check that the total equals occupied + empty
    EXPECT_EQ(grid.getTotalVoxelCount(), 
              grid.getOccupiedVoxelCount() + grid.getEmptyVoxelCount());
}

TEST_F(VoxelGridTest, MultipleAtoms) {
    // Create multiple atoms
    std::vector<std::unique_ptr<Atom>> atoms;
    
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    atoms.push_back(std::make_unique<Atom>("N", 0.56, 14.007));
    atoms.back()->setCoordinates(2.0, 0.0, 0.0);
    atoms.back()->setId(1);
    
    atoms.push_back(std::make_unique<Atom>("O", 0.48, 15.999));
    atoms.back()->setCoordinates(0.0, 2.0, 0.0);
    atoms.back()->setId(2);
    
    // Create grid
    VoxelGrid grid(atoms, 0.5, 1.0);
    
    // Check that multiple voxels are occupied
    EXPECT_GT(grid.getOccupiedVoxelCount(), 2);
}

TEST_F(VoxelGridTest, DifferentVoxelSizes) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    // Test different voxel sizes
    VoxelGrid grid1(atoms, 0.5, 1.0);
    VoxelGrid grid2(atoms, 1.0, 1.0);
    VoxelGrid grid3(atoms, 2.0, 1.0);
    
    // Larger voxel size should result in fewer total voxels
    EXPECT_GT(grid1.getTotalVoxelCount(), grid2.getTotalVoxelCount());
    EXPECT_GT(grid2.getTotalVoxelCount(), grid3.getTotalVoxelCount());
}

TEST_F(VoxelGridTest, InvalidVoxelSizeThrows) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    // Negative voxel size should throw
    EXPECT_THROW(VoxelGrid(atoms, -1.0, 0.0), std::invalid_argument);
    
    // Zero voxel size should throw
    EXPECT_THROW(VoxelGrid(atoms, 0.0, 0.0), std::invalid_argument);
}

TEST_F(VoxelGridTest, GridDimensions) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 1.0, 1.0);
    
    const auto& dims = grid.getDimensions();
    
    // All dimensions should be positive
    EXPECT_GT(dims[0], 0);
    EXPECT_GT(dims[1], 0);
    EXPECT_GT(dims[2], 0);
    
    // Total voxels should equal product of dimensions
    EXPECT_EQ(grid.getTotalVoxelCount(), dims[0] * dims[1] * dims[2]);
}

TEST_F(VoxelGridTest, OccupiedVoxelProperties) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(5.0, 5.0, 5.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 1.0, 1.0);
    
    const auto& occupiedVoxels = grid.getOccupiedVoxels();
    
    // Check that occupied voxels have valid properties
    for (const auto& voxel : occupiedVoxels) {
        EXPECT_TRUE(voxel.occupied);
        EXPECT_GT(voxel.atomIds.size(), 0);
        
        // Check that voxel bounds are consistent
        EXPECT_LT(voxel.min.x, voxel.max.x);
        EXPECT_LT(voxel.min.y, voxel.max.y);
        EXPECT_LT(voxel.min.z, voxel.max.z);
        
        // Check that center is actually in the middle
        double expectedCenterX = (voxel.min.x + voxel.max.x) * 0.5;
        double expectedCenterY = (voxel.min.y + voxel.max.y) * 0.5;
        double expectedCenterZ = (voxel.min.z + voxel.max.z) * 0.5;
        
        EXPECT_NEAR(voxel.center.x, expectedCenterX, 1e-10);
        EXPECT_NEAR(voxel.center.y, expectedCenterY, 1e-10);
        EXPECT_NEAR(voxel.center.z, expectedCenterZ, 1e-10);
    }
}

// VoxelMeshGenerator tests
TEST_F(VoxelMeshGeneratorTest, BasicMeshGeneration) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 1.0, 1.0);
    HexMesh mesh = VoxelMeshGenerator::generateHexMesh(grid);
    
    // Mesh should have nodes and elements
    EXPECT_GT(mesh.getNodeCount(), 0);
    EXPECT_GT(mesh.getElementCount(), 0);
    
    // Number of elements should equal occupied voxels
    EXPECT_EQ(mesh.getElementCount(), grid.getOccupiedVoxelCount());
}

TEST_F(VoxelMeshGeneratorTest, NodeDeduplication) {
    // Create two atoms close to each other to create adjacent voxels
    std::vector<std::unique_ptr<Atom>> atoms;
    
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(5.0, 5.0, 5.0);
    atoms.back()->setId(0);
    
    atoms.push_back(std::make_unique<Atom>("N", 0.56, 14.007));
    atoms.back()->setCoordinates(5.5, 5.0, 5.0);
    atoms.back()->setId(1);
    
    VoxelGrid grid(atoms, 0.5, 1.0);
    HexMesh mesh = VoxelMeshGenerator::generateHexMesh(grid);
    
    // Node count should be less than 8 * element count due to sharing
    EXPECT_LT(mesh.getNodeCount(), mesh.getElementCount() * 8);
}

TEST_F(VoxelMeshGeneratorTest, EmptyGrid) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.01, 12.011)); // Very small radius
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 10.0, 1.0); // Very large voxels
    HexMesh mesh = VoxelMeshGenerator::generateHexMesh(grid);
    
    // If no voxels are occupied, mesh should be empty
    if (grid.getOccupiedVoxelCount() == 0) {
        EXPECT_EQ(mesh.getNodeCount(), 0);
        EXPECT_EQ(mesh.getElementCount(), 0);
    }
}

TEST_F(VoxelMeshGeneratorTest, ElementConnectivity) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(0.0, 0.0, 0.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 1.0, 1.0);
    HexMesh mesh = VoxelMeshGenerator::generateHexMesh(grid);
    
    // Check that all element connectivity indices are valid
    for (const auto& element : mesh.elements) {
        for (int nodeIdx : element) {
            EXPECT_GE(nodeIdx, 0);
            EXPECT_LT(nodeIdx, static_cast<int>(mesh.getNodeCount()));
        }
    }
}

TEST_F(VoxelMeshGeneratorTest, StandardHexOrdering) {
    std::vector<std::unique_ptr<Atom>> atoms;
    atoms.push_back(std::make_unique<Atom>("C", 0.67, 12.011));
    atoms.back()->setCoordinates(5.0, 5.0, 5.0);
    atoms.back()->setId(0);
    
    VoxelGrid grid(atoms, 0.5, 1.0);  // Use smaller voxels to ensure occupancy
    HexMesh mesh = VoxelMeshGenerator::generateHexMesh(grid);
    
    // Should have at least one element
    ASSERT_GT(mesh.getElementCount(), 0);
    
    // Check first element has 8 nodes
    const auto& element = mesh.elements[0];
    EXPECT_EQ(element.size(), 8);
    
    // All 8 node indices should be valid
    for (size_t i = 0; i < 8; ++i) {
        EXPECT_GE(element[i], 0);
        EXPECT_LT(element[i], static_cast<int>(mesh.getNodeCount()));
    }
}

TEST_F(VoxelMeshGeneratorTest, ParallelConsistency) {
    // Test that parallel execution produces consistent results
    // Create a moderately sized structure for meaningful parallelization
    std::vector<std::unique_ptr<Atom>> atoms;
    
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            auto atom = std::make_unique<Atom>("C", 1.0, 12.011);
            atom->setCoordinates(i * 2.0, j * 2.0, 5.0);
            atom->setId(i * 10 + j);
            atoms.push_back(std::move(atom));
        }
    }
    
    VoxelGrid grid(atoms, 0.8, 1.0);
    
    // Generate mesh multiple times - should be consistent
    HexMesh mesh1 = VoxelMeshGenerator::generateHexMesh(grid);
    HexMesh mesh2 = VoxelMeshGenerator::generateHexMesh(grid);
    
    // Verify consistency
    EXPECT_EQ(mesh1.getNodeCount(), mesh2.getNodeCount());
    EXPECT_EQ(mesh1.getElementCount(), mesh2.getElementCount());
    
    // Verify all elements are identical
    ASSERT_EQ(mesh1.elements.size(), mesh2.elements.size());
    for (size_t i = 0; i < mesh1.elements.size(); ++i) {
        for (size_t j = 0; j < 8; ++j) {
            EXPECT_EQ(mesh1.elements[i][j], mesh2.elements[i][j]);
        }
    }
}

