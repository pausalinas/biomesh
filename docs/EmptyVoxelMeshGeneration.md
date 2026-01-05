# Empty Voxel Mesh Generation with GiD Export

This module provides functionality to generate hexahedral finite element meshes from **empty voxels** (non-occupied regions) and export them to GiD format. This is complementary to the standard `VoxelMeshGenerator` which processes occupied voxels.

## Overview

The `EmptyVoxelMeshGenerator` creates hexahedral meshes from the void space around molecular structures. This is particularly useful for:

- **Computational Fluid Dynamics (CFD)**: Modeling fluid flow around proteins or molecules
- **Structural Analysis**: Analyzing the external environment of molecular structures
- **Void Space Analysis**: Studying empty regions in molecular assemblies
- **Inverse Modeling**: Creating negative space representations

## Use Cases

### 1. CFD Simulations
Empty voxel meshes represent the fluid domain surrounding a molecule, enabling:
- Flow field analysis around proteins
- Drag and lift coefficient calculations
- Mass transfer studies

### 2. External Structural Analysis
Model the surrounding medium for:
- Pressure distribution analysis
- Thermal analysis of the environment
- Mechanical interaction with external forces

### 3. Void Space Characterization
Quantify and visualize:
- Cavity volumes in protein assemblies
- Pore networks in molecular frameworks
- Accessible surface areas

## GiD File Format

GiD is a universal pre/post-processing software widely used in finite element and computational fluid dynamics applications. The `.msh` file format follows this structure:

```
MESH dimension 3 ElemType Hexahedra Nnode 8

Coordinates
1 x1 y1 z1
2 x2 y2 z2
...
N xN yN zN
End Coordinates

Elements
1 n1 n2 n3 n4 n5 n6 n7 n8
2 n1 n2 n3 n4 n5 n6 n7 n8
...
M n1 n2 n3 n4 n5 n6 n7 n8
End Elements
```

**Key Features:**
- **1-based indexing**: Node and element indices start at 1 (not 0)
- **Fixed precision**: Coordinates use 6 decimal places
- **Standard ordering**: Hexahedral elements follow standard node ordering

## Usage

### Basic Usage

```cpp
#include "biomesh/EmptyVoxelMeshGenerator.hpp"
#include "biomesh/GiDExporter.hpp"
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/PDBParser.hpp"
#include "biomesh/AtomBuilder.hpp"

// Parse and enrich atoms
auto basicAtoms = PDBParser::parsePDBFile("protein.pdb");
AtomBuilder builder;
auto enrichedAtoms = builder.buildAtoms(basicAtoms);

// Create voxel grid
double voxelSize = 1.0;  // 1.0 Angstrom voxels
double padding = 2.0;     // 2.0 Angstrom padding
VoxelGrid voxelGrid(enrichedAtoms, voxelSize, padding);

// Generate hexahedral mesh from empty voxels
HexMesh mesh = EmptyVoxelMeshGenerator::generateHexMesh(voxelGrid);

// Export to GiD format
bool success = GiDExporter::exportToGiD(mesh, "empty_mesh.msh");
```

### Command-Line Tool

A dedicated command-line tool is provided for easy mesh generation:

```bash
# Basic usage
./empty_voxel_to_gid protein.pdb 1.0 empty_mesh.msh

# With custom padding
./empty_voxel_to_gid protein.pdb 1.0 empty_mesh.msh 2.5

# Larger voxels for faster processing
./empty_voxel_to_gid protein.pdb 2.0 empty_mesh.msh 3.0
```

**Arguments:**
- `pdb_file`: Path to input PDB file
- `voxel_size`: Edge length of voxels in Angstroms
- `output_file`: Output GiD mesh file path (`.msh` extension)
- `padding`: Optional padding around bounding box (default: 2.0 Å)

### Example Output

```
BioMesh - Empty Voxel Mesh Generator with GiD Export
======================================================

Loading PDB file: protein.pdb
  Loaded 245 atoms

Enriching atoms with physical properties...
  Enriched 245 atoms

Creating voxel grid...
  Voxel size: 1.0 Å
  Padding: 2.0 Å

=== Voxel Grid Statistics ===
Voxel size: 1.0 Å
Grid dimensions: 45 x 52 x 38
Total voxels: 88920
Occupied voxels: 1234 (1.4%)
Empty voxels: 87686 (98.6%)
=============================

Generating hexahedral mesh from empty voxels...
  Generated mesh:
    Nodes: 245678
    Elements: 87686
    Node sharing efficiency: 65.2%

Exporting to GiD format: empty_mesh.msh
  Export successful!

Mesh file written to: empty_mesh.msh
```

## Node Ordering Convention

Both `EmptyVoxelMeshGenerator` and `VoxelMeshGenerator` use the same standard hexahedral element node ordering:

```
       7--------6
      /|       /|
     / |      / |
    4--------5  |
    |  3-----|--2
    | /      | /
    |/       |/
    0--------1
```

- Node 0: (min.x, min.y, min.z) - bottom-left-back
- Node 1: (max.x, min.y, min.z) - bottom-right-back
- Node 2: (max.x, max.y, min.z) - bottom-right-front
- Node 3: (min.x, max.y, min.z) - bottom-left-front
- Node 4: (min.x, min.y, max.z) - top-left-back
- Node 5: (max.x, min.y, max.z) - top-right-back
- Node 6: (max.x, max.y, max.z) - top-right-front
- Node 7: (min.x, max.y, max.z) - top-left-front

This ordering is compatible with most FEM software including GiD, ANSYS, and Abaqus.

## Compilation

The module is automatically built with BioMesh:

```bash
mkdir build && cd build
cmake ..
make -j4
```

The example executable will be created as `build/empty_voxel_to_gid`.

## Performance Considerations

### Mesh Size Warning

**Important**: Empty voxel meshes are typically **much larger** than occupied voxel meshes:
- Occupied voxels: typically 1-10% of total voxels
- Empty voxels: typically 90-99% of total voxels

For a grid with 100,000 voxels:
- Occupied mesh: ~5,000 elements
- Empty mesh: ~95,000 elements (19x larger!)

### Recommendations

1. **Start with larger voxel sizes** (2.0-3.0 Å) for initial testing
2. **Use smaller voxel sizes** (0.5-1.0 Å) only when high resolution is needed
3. **Monitor file sizes**: Large meshes can produce multi-gigabyte files
4. **Consider memory**: Node deduplication uses hash maps that scale with mesh size

### OpenMP Parallelization

The mesh generator uses OpenMP for parallel corner node computation:
- Enabled automatically if OpenMP is available during compilation
- Provides significant speedup for large meshes
- Scales well with number of CPU cores

### Voxel Size Guidelines

| Voxel Size | Use Case | Typical Elements | Processing Time |
|------------|----------|------------------|-----------------|
| 3.0 Å | Quick preview | < 10,000 | Seconds |
| 2.0 Å | Standard CFD | 10,000 - 50,000 | Seconds |
| 1.0 Å | High resolution | 50,000 - 500,000 | Minutes |
| 0.5 Å | Very fine detail | 500,000+ | Hours |

## Error Handling

The implementation handles several error conditions:

### Empty Voxel Grid
```cpp
HexMesh mesh = EmptyVoxelMeshGenerator::generateHexMesh(voxelGrid);
if (mesh.getElementCount() == 0) {
    // No empty voxels found
}
```

### File I/O Errors
```cpp
bool success = GiDExporter::exportToGiD(mesh, "output.msh");
if (!success) {
    // File could not be written
}
```

### Invalid Parameters
The command-line tool validates:
- Positive voxel size
- Non-negative padding
- Valid input file path
- Writable output directory

## Integration with Existing Code

The module seamlessly integrates with existing BioMesh functionality:

```cpp
// Use the same VoxelGrid for both occupied and empty meshes
VoxelGrid voxelGrid(enrichedAtoms, voxelSize, padding);

// Generate occupied mesh (molecules)
HexMesh occupiedMesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);

// Generate empty mesh (surrounding space)
HexMesh emptyMesh = EmptyVoxelMeshGenerator::generateHexMesh(voxelGrid);

// Export both meshes
GiDExporter::exportToGiD(occupiedMesh, "molecule.msh");
GiDExporter::exportToGiD(emptyMesh, "fluid_domain.msh");
```

## Advanced Examples

### Multi-Resolution Analysis
```cpp
std::vector<double> resolutions = {3.0, 2.0, 1.0};
for (double voxelSize : resolutions) {
    VoxelGrid grid(atoms, voxelSize, 2.0);
    HexMesh mesh = EmptyVoxelMeshGenerator::generateHexMesh(grid);
    
    std::string filename = "mesh_" + std::to_string(voxelSize) + ".msh";
    GiDExporter::exportToGiD(mesh, filename);
}
```

### Combined Occupied and Empty Meshes
```cpp
// Create grids with different padding for interior and exterior
VoxelGrid interiorGrid(atoms, 1.0, 0.5);  // Tight fit for molecule
VoxelGrid exteriorGrid(atoms, 1.0, 5.0);  // Large domain for CFD

HexMesh molecule = VoxelMeshGenerator::generateHexMesh(interiorGrid);
HexMesh fluid = EmptyVoxelMeshGenerator::generateHexMesh(exteriorGrid);

GiDExporter::exportToGiD(molecule, "molecule.msh");
GiDExporter::exportToGiD(fluid, "fluid_domain.msh");
```

## References

- GiD Software: https://www.gidsimulation.com/
- VoxelGrid documentation: See `docs/VoxelizationMeshGeneration.md`
- HexMesh data structure: See `include/biomesh/HexMesh.hpp`

## See Also

- `VoxelMeshGenerator`: For generating meshes from occupied voxels
- `VoxelGrid`: For voxel grid creation and management
- `GiDExporter`: For exporting meshes to GiD format
