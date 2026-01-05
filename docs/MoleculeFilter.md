# MoleculeFilter Documentation

## Overview

The MoleculeFilter component provides a flexible way to filter PDB files by biomolecule type before processing. This is useful when you want to focus on specific molecule types (e.g., proteins only) while excluding others (e.g., water molecules, ions).

## Features

- **Residue Classification**: Automatically identifies proteins, nucleic acids (DNA/RNA), water, ions, and other molecules
- **Preset Filters**: Common filtering scenarios ready to use
- **Custom Filters**: Builder pattern for creating custom filter combinations
- **Preserves Atom Properties**: Filtered atoms retain all their properties (coordinates, element, residue info, etc.)

## Quick Start

### Basic Usage

```cpp
#include "biomesh/PDBParser.hpp"
#include "biomesh/MoleculeFilter.hpp"

using namespace biomesh;

// Parse PDB file
auto atoms = PDBParser::parsePDBFile("complex.pdb");

// Filter to keep only protein atoms
auto filter = MoleculeFilter::proteinOnly();
auto proteinAtoms = filter.filter(atoms);
```

### Preset Filters

The library provides several preset filters for common scenarios:

```cpp
// Keep everything (no filtering)
auto filterAll = MoleculeFilter::all();

// Keep only protein atoms
auto filterProtein = MoleculeFilter::proteinOnly();

// Keep only nucleic acids (DNA/RNA)
auto filterNucleic = MoleculeFilter::nucleicAcidOnly();

// Remove water molecules
auto filterNoWater = MoleculeFilter::noWater();
```

### Custom Filters

Use the builder pattern to create custom filter combinations:

```cpp
MoleculeFilter customFilter;
customFilter.setKeepProteins(true)
            .setKeepNucleicAcids(true)
            .setKeepWater(false)
            .setKeepIons(false)
            .setKeepOthers(false);

auto filtered = customFilter.filter(atoms);
```

## ResidueClassifier

The `ResidueClassifier` class provides static methods to classify residues:

```cpp
#include "biomesh/ResidueClassifier.hpp"

// Check if residue is a protein
bool isProtein = ResidueClassifier::isProtein("ALA");  // true

// Check if residue is DNA
bool isDNA = ResidueClassifier::isDNA("DA");  // true

// Check if residue is RNA
bool isRNA = ResidueClassifier::isRNA("A");  // true

// Check if residue is nucleic acid (DNA or RNA)
bool isNucleic = ResidueClassifier::isNucleicAcid("DG");  // true

// Check if residue is water
bool isWater = ResidueClassifier::isWater("HOH");  // true

// Check if residue is an ion
bool isIon = ResidueClassifier::isIon("CA");  // true
```

### Supported Residue Types

**Proteins (Standard Amino Acids):**
- ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

**Proteins (Non-Standard Amino Acids):**
- MSE (Selenomethionine)
- SEC (Selenocysteine)
- PYL (Pyrrolysine)

**DNA:**
- DA (Deoxyadenosine)
- DT (Deoxythymidine)
- DG (Deoxyguanosine)
- DC (Deoxycytidine)
- DU (Deoxyuridine)

**RNA:**
- A (Adenosine)
- U (Uridine)
- G (Guanosine)
- C (Cytidine)
- ADE, URA, GUA, CYT (Alternative notations)

**Water:**
- HOH, WAT, H2O, SOL, TIP, TIP3, TIP4

**Common Ions:**
- NA (Sodium)
- CL (Chloride)
- K (Potassium)
- CA (Calcium)
- MG (Magnesium)
- ZN (Zinc)
- FE (Iron)
- CU (Copper)
- MN (Manganese)

## API Reference

### MoleculeFilter Class

#### Static Factory Methods

```cpp
static MoleculeFilter all();
static MoleculeFilter proteinOnly();
static MoleculeFilter nucleicAcidOnly();
static MoleculeFilter noWater();
```

#### Builder Methods

```cpp
MoleculeFilter& setKeepProteins(bool keep);
MoleculeFilter& setKeepNucleicAcids(bool keep);
MoleculeFilter& setKeepWater(bool keep);
MoleculeFilter& setKeepIons(bool keep);
MoleculeFilter& setKeepOthers(bool keep);
```

#### Filtering Methods

```cpp
// Filter a collection of atoms
std::vector<std::unique_ptr<Atom>> filter(const std::vector<std::unique_ptr<Atom>>& atoms) const;

// Check if a single atom should be kept
bool shouldKeep(const Atom& atom) const;
```

### ResidueClassifier Class

All methods are static:

```cpp
static bool isProtein(const std::string& residueName);
static bool isNucleicAcid(const std::string& residueName);
static bool isDNA(const std::string& residueName);
static bool isRNA(const std::string& residueName);
static bool isWater(const std::string& residueName);
static bool isIon(const std::string& residueName);
```

## Complete Example

```cpp
#include "biomesh/BioMesh.hpp"
#include <iostream>

using namespace biomesh;

int main() {
    try {
        // Parse PDB file (includes ATOM and HETATM records)
        auto atoms = PDBParser::parsePDBFile("1abc.pdb");
        
        std::cout << "Total atoms: " << atoms.size() << std::endl;
        
        // Create custom filter: proteins and DNA only, no water or ions
        MoleculeFilter filter;
        filter.setKeepProteins(true)
              .setKeepNucleicAcids(true)  // Includes DNA and RNA
              .setKeepWater(false)
              .setKeepIons(false)
              .setKeepOthers(false);
        
        auto filtered = filter.filter(atoms);
        
        std::cout << "Filtered atoms: " << filtered.size() << std::endl;
        
        // Enrich with physical properties
        AtomBuilder builder;
        auto enriched = builder.buildAtoms(filtered);
        
        // Calculate bounding box
        BoundingBox bbox(enriched, 2.0);
        
        std::cout << "Volume: " << bbox.getVolume() << " Ų" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
```

## Integration with Workflow

The MoleculeFilter fits into the BioMesh pipeline as an optional preprocessing step:

```
PDB File → PDBParser → [MoleculeFilter] → AtomBuilder → VoxelGrid → HexMesh
                            ↑
                      (optional preprocessing)
```

Typical workflow:

1. **Parse PDB**: Extract all atoms from the file
2. **Filter (Optional)**: Select desired molecule types
3. **Enrich**: Add physical properties (radius, mass)
4. **Voxelize**: Create voxel grid
5. **Generate Mesh**: Create hexahedral mesh

## Notes

- Classification is case-insensitive (e.g., "ALA", "ala", and "Ala" are all recognized)
- Unknown residues (not in any predefined category) are kept when `setKeepOthers(true)`
- Filtering creates deep copies of atoms, preserving all properties
- The parser now handles both ATOM and HETATM records from PDB files

## Testing

Run the included demo to test filtering:

```bash
cd build
./filter_demo                          # Uses default test file
./filter_demo ../data/mixed_molecules.pdb  # Test with mixed molecules
```

The test suite includes comprehensive tests for:
- ResidueClassifier: All molecule type classifications
- MoleculeFilter: Preset filters, custom filters, edge cases
- Atom property preservation during filtering
