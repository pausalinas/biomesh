# BioMesh - Molecular Mesh Generator

[![License:  MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![C++17](https://img.shields.io/badge/C++-17-blue.svg)](https://isocpp.org/)
[![CMake](https://img.shields.io/badge/CMake-3.14+-064F8C. svg)](https://cmake.org/)

A modern C++ toolkit for generating hexahedral finite element meshes from molecular structures in PDB format. BioMesh2 provides both **occupied voxel meshes** (molecule) and **empty voxel meshes** (void space) for CFD, FEM, and molecular dynamics applications.

---

## Features

### Core Capabilities
- **PDB Parsing** - Read and extract atom coordinates from Protein Data Bank files
- **Atomic Enrichment** - Automatically assign van der Waals radii and atomic masses
- **Bounding Box Calculation** - Compute 3D molecular bounding boxes with padding
- **Uniform Voxelization** - Tessellate space into regular cubic voxels
- **Dual Mesh Generation** - Create meshes from both occupied and empty voxels
- **GiD Export** - Export to industry-standard GiD mesh format (. msh)
- **Molecule Filtering** - Filter by protein, nucleic acids, water, ions, etc.

### Advanced Features
- **OpenMP Parallelization** - Multi-threaded mesh generation for large structures
- **Command-Line Tool** - Easy-to-use bash script with config file support
- **High Precision** - Configurable voxel resolution from 3. 0 Å to 0.01 Å
- **Comprehensive Testing** - 60+ unit tests covering all functionality
- **Modern C++17** - RAII, smart pointers, STL containers

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/pausalinas/BioMesh2.git
cd BioMesh2

# Build
mkdir build && cd build
cmake .. 
make -j4

# Verify
./occupied_voxel_to_gid
./empty_voxel_to_gid
```

### Basic Usage

```bash
# Return to project root
cd ..

# Generate both occupied and empty meshes
./biomesh2. sh protein.pdb

# Output:
#   mesh_occupied.msh  (molecule mesh)
#   mesh_empty.msh     (fluid domain mesh)
```

### High-Resolution Example

```bash
# Download a real protein
wget https://files.rcsb.org/download/1CRN.pdb

# Generate high-resolution meshes
./biomesh2.sh 1CRN.pdb -v 0.5 -o crambin_highres --verbose

# Outputs:
#   crambin_highres_occupied.msh
#   crambin_highres_empty.msh
```

---

## Documentation

- **[Getting Started Guide](docs/getting_started.md)** - Tutorial for new users
- **[API Reference](docs/api_reference.md)** - Complete C++ API documentation
- **[Configuration Guide](docs/configuration.md)** - Config file reference
- **[Examples](docs/examples.md)** - Practical use cases and workflows
- **[Voxelization Details](docs/VoxelizationMeshGeneration.md)** - Technical documentation
- **[Empty Mesh Generation](docs/EmptyVoxelMeshGeneration.md)** - CFD-specific guide

---

## Project Structure

```
BioMesh2/
├── include/biomesh/              # Public headers
│   ├── Atom.hpp                  # Atom class with physical properties
│   ├── AtomicSpec.hpp            # Atomic database (radii, masses)
│   ├── PDBParser.hpp             # PDB file parser
│   ├── AtomBuilder.hpp           # Atom enrichment
│   ├── BoundingBox.hpp           # 3D bounding box calculator
│   ├── VoxelGrid.hpp             # Voxel grid data structure
│   ├── VoxelMeshGenerator.hpp    # Occupied voxel → mesh
│   ├── EmptyVoxelMeshGenerator.hpp # Empty voxel → mesh
│   ├── GiDExporter.hpp           # GiD format exporter
│   ├── MoleculeFilter.hpp        # Molecule type filtering
│   ├── ResidueClassifier.hpp     # Residue classification
│   └── BioMesh.hpp               # Main convenience header
├── src/                          # Implementation files
├── examples/                     # Example programs
│   ├── occupied_voxel_to_gid.cpp # Occupied mesh generator
│   ├── empty_voxel_to_gid.cpp    # Empty mesh generator
│   ├── main.cpp                  # Basic demo
│   ├── voxel_demo.cpp            # Voxelization demo
│   ├── filter_demo.cpp           # Filtering demo
│   └── filter_workflow_example.cpp
├── tests/                        # Unit tests (60+ tests)
├── data/                         # Sample PDB files
├── docs/                         # Documentation
├── config/                       # Configuration files
├── biomesh2.sh                   # Command-line tool
└── CMakeLists.txt                # Build configuration
```

---

## Use Cases

### Computational Fluid Dynamics (CFD)
Generate fluid domain meshes around proteins for flow analysis: 
```bash
./biomesh2.sh protein.pdb --occupied false -v 0.5 -p 10. 0 -o fluid_domain
```

### Finite Element Analysis (FEM)
Create molecular meshes for structural analysis:
```bash
./biomesh2.sh protein.pdb --empty false -v 1.0 -o structural_mesh
```

### Multi-Physics Simulations
Generate both molecule and surrounding medium: 
```bash
./biomesh2.sh protein.pdb -v 0.5 -o simulation
# Outputs both occupied and empty meshes
```

### Convergence Studies
Multi-resolution mesh generation: 
```bash
for voxel in 2. 0 1.0 0.5 0.25; do
    ./biomesh2.sh protein.pdb -v $voxel -o mesh_${voxel}
done
```

---

## Command-Line Tool

### Basic Syntax

```bash
./biomesh2.sh [OPTIONS] <input. pdb>
```

### Key Options

| Option | Description | Default |
|--------|-------------|---------|
| `-v, --voxel-size VALUE` | Voxel edge length (Å) | 1.0 |
| `-p, --padding VALUE` | Bounding box padding (Å) | 2.0 |
| `-o, --output BASENAME` | Output file basename | mesh |
| `--occupied BOOL` | Generate occupied mesh | true |
| `--empty BOOL` | Generate empty mesh | true |
| `--filter TYPE` | Molecule filter | none |
| `-c, --config FILE` | Use configuration file | - |
| `--batch LIST` | Process multiple files | - |
| `-V, --verbose` | Verbose output | false |

### Configuration File

```bash
# Generate template
./biomesh2.sh --generate-config > my_config.conf

# Edit configuration
nano my_config.conf

# Run with config
./biomesh2.sh --config my_config.conf
```

**Example config:**
```ini
[input]
file = protein.pdb

[voxelization]
voxel_size = 0.5
padding = 2.0

[output]
basename = high_res_mesh
generate_occupied = true
generate_empty = true

[options]
verbose = true
```

---

## C++ API

### Complete Workflow

```cpp
#include "biomesh/BioMesh.hpp"
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/VoxelMeshGenerator.hpp"
#include "biomesh/EmptyVoxelMeshGenerator.hpp"
#include "biomesh/GiDExporter.hpp"

using namespace biomesh;

int main() {
    // Parse PDB
    auto basicAtoms = PDBParser::parsePDBFile("protein.pdb");
    
    // Enrich with physical properties
    AtomBuilder builder;
    auto enrichedAtoms = builder.buildAtoms(basicAtoms);
    
    // Create voxel grid
    VoxelGrid voxelGrid(enrichedAtoms, 1.0, 2.0); // 1.0 Å voxels, 2.0 Å padding
    
    // Generate occupied mesh (molecule)
    HexMesh occupiedMesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);
    
    // Generate empty mesh (void space)
    HexMesh emptyMesh = EmptyVoxelMeshGenerator:: generateHexMesh(voxelGrid);
    
    // Export to GiD format
    GiDExporter:: exportToGiD(occupiedMesh, "occupied. msh");
    GiDExporter::exportToGiD(emptyMesh, "empty.msh");
    
    std::cout << "Occupied:  " << occupiedMesh.getElementCount() << " elements\n";
    std::cout << "Empty: " << emptyMesh. getElementCount() << " elements\n";
    
    return 0;
}
```

### Molecule Filtering

```cpp
#include "biomesh/MoleculeFilter.hpp"

// Filter to protein only
auto filter = MoleculeFilter::proteinOnly();
auto proteinAtoms = filter.filter(atoms);

// Custom filter
MoleculeFilter customFilter;
customFilter.setKeepProteins(true)
            .setKeepWater(false)
            .setKeepIons(false);
auto filtered = customFilter.filter(atoms);
```

---

## Resolution Guidelines

| Voxel Size | Resolution | Elements (typical) | Use Case |
|------------|-----------|-------------------|----------|
| 3.0 Å | Very Low | 1K - 10K | Quick preview |
| 2.0 Å | Low | 10K - 50K | Draft visualization |
| 1.0 Å | Medium | 50K - 200K | Standard work |
| 0.5 Å | High | 200K - 1M | Production quality |
| 0.25 Å | Very High | 1M - 10M | Detailed analysis |
| 0.1 Å | Ultra High | 10M+ | Research/HPC |

**Performance Tip:** Start with 2.0 Å for testing, then refine to 0.5-1.0 Å for production. 

---

## Build Options

### Prerequisites

- **C++17 compiler** (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake 3.14+**
- **GoogleTest** (optional, for tests)
- **OpenMP** (optional, for parallelization)

### Build Configuration

```bash
# Standard build
cmake .. 
make -j4

# Release build (optimized)
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4

# With OpenMP parallelization
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-O3 -fopenmp" .. 
make -j4

# Run tests
make test
```

### OpenMP Parallelization

```bash
# Set thread count
export OMP_NUM_THREADS=16

# Run with parallelization
./occupied_voxel_to_gid protein.pdb 0.5 mesh. msh 2. 0
```

---

## Testing

Comprehensive test suite with 60+ tests: 

```bash
cd build
./biomesh_tests

# Or via CMake
make test
```

**Test Coverage:**
- Atom construction and properties
- Atomic database operations
- PDB parsing (ATOM and HETATM records)
- Bounding box calculations
- Voxel grid generation
- Mesh generation (occupied & empty)
- Molecule filtering
- GiD export format

---

## Output Formats

### GiD Mesh Format (. msh)

Industry-standard format for FEM/CFD software:

```
MESH dimension 3 ElemType Hexahedra Nnode 8

Coordinates
1 x1 y1 z1
2 x2 y2 z2
... 
End Coordinates

Elements
1 n1 n2 n3 n4 n5 n6 n7 n8
2 n1 n2 n3 n4 n5 n6 n7 n8
... 
End Elements
```

**Compatible Software:**
- GiD (pre/post-processing)
- ANSYS
- Abaqus
- CalculiX
- Custom FEM/CFD codes

---

## 🎓 Examples

### Example 1: Water Molecule
```bash
cat > water. pdb <<EOF
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
ATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00  0.00           H
ATOM      3  H2  HOH A   1      -0.240   0.927   0.000  1.00  0.00           H
END
EOF

./biomesh2.sh water. pdb -v 0.3 -p 1.0 -o water_mesh
```

### Example 2: Batch Processing
```bash
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/1MSO.pdb
wget https://files.rcsb.org/download/2LYZ.pdb

./biomesh2.sh --batch 1CRN.pdb,1MSO.pdb,2LYZ.pdb -v 1.0
```

### Example 3: CFD Preparation
```bash
./biomesh2.sh protein.pdb \
  --occupied false \
  -v 0.5 \
  -p 15.0 \
  -o cfd_domain \
  --verbose
```

---

## Supported Elements

Built-in atomic database with van der Waals radii and masses:

| Element | Radius (Å) | Mass (Da) | Category |
|---------|-----------|-----------|----------|
| H | 1.20 | 1.008 | Main |
| C | 1.70 | 12.011 | Main |
| N | 1.55 | 14.007 | Main |
| O | 1.52 | 15.999 | Main |
| P | 1.80 | 30.974 | Main |
| S | 1.80 | 32.065 | Main |
| F, Cl, Br, I | varies | varies | Halogens |
| Na, Mg, K, Ca | varies | varies | Metals |
| Fe, Zn, Se | varies | varies | Trace elements |

Additional elements can be easily added to the database.

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

**Development Guidelines:**
- Follow C++17 best practices
- Add unit tests for new features
- Update documentation
- Ensure all tests pass

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- Protein Data Bank (PDB) for molecular structure format
- GiD development team for mesh format specification
- OpenMP for parallelization support
- GoogleTest for testing framework

---

## 📞 Support

- **Documentation:** [docs/](docs/)
- **Issues:** [GitHub Issues](https://github.com/pausalinas/BioMesh/issues)
- **Discussions:** [GitHub Discussions](https://github.com/pausalinas/BioMesh/discussions)

---

## Roadmap

- [ ] VTK export format
- [ ] Python bindings (pybind11)
- [ ] Adaptive mesh refinement
- [ ] Surface mesh extraction
- [ ] Tetrahedral mesh generation
- [ ] ParaView plugin

---

## Citation

If you use BioMesh in your research, please cite:

```bibtex
@software{biomesh,
  author = {Salinas, Paulina},
  title = {BioMesh: Molecular Hexahedral Mesh Generator},
  year = {2026},
  url = {https://github.com/pausalinas/BioMesh}
}
```

---

**Made with ᡣ𐭩 for the computational biology and CFD communities**
