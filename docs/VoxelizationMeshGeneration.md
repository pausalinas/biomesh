# Voxelization-Based Hexahedral Mesh Generation

This module provides functionality to generate hexahedral finite element meshes using a uniform voxelization approach. The implementation tessellates the molecular bounding box into uniform cubic voxels and creates hexahedral elements only for voxels that intersect with atom spheres.

## Features

- **Uniform Voxelization**: Tessellates the entire bounding box into regular cubic voxels
- **User-Configurable Voxel Size**: Allows precise control over mesh resolution via voxel edge length
- **Occupancy Detection**: Classifies voxels based on intersection with atom spheres (center-point method)
- **Separate Storage**: Maintains both occupied and empty voxels for potential future use
- **Node Deduplication**: Automatically detects and reuses shared corner nodes between adjacent voxels
- **Standard Element Ordering**: Uses standard hexahedral element node ordering compatible with FEM software
- **Memory Efficient**: Uses hash maps for O(1) node lookup and deduplication

## Usage

### Basic Voxelization

```cpp
#include "biomesh/VoxelGrid.hpp"
#include "biomesh/VoxelMeshGenerator.hpp"
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

// Generate hexahedral mesh
HexMesh mesh = VoxelMeshGenerator::generateHexMesh(voxelGrid);

std::cout << "Nodes: " << mesh.getNodeCount() << "\n";
std::cout << "Elements: " << mesh.getElementCount() << "\n";
```

### Using Existing Bounding Box

```cpp
// Compute bounding box separately
BoundingBox bbox(enrichedAtoms, padding);

// Create voxel grid with existing bounding box
VoxelGrid voxelGrid(bbox, enrichedAtoms, voxelSize);
```

### Accessing Grid Information

```cpp
// Get grid dimensions
const auto& dims = voxelGrid.getDimensions();
std::cout << "Grid: " << dims[0] << " x " << dims[1] << " x " << dims[2] << "\n";

// Get voxel counts
std::cout << "Total voxels: " << voxelGrid.getTotalVoxelCount() << "\n";
std::cout << "Occupied voxels: " << voxelGrid.getOccupiedVoxelCount() << "\n";
std::cout << "Empty voxels: " << voxelGrid.getEmptyVoxelCount() << "\n";

// Print detailed statistics
voxelGrid.printStatistics();
```

### Working with Individual Voxels

```cpp
// Access occupied voxels
const auto& occupiedVoxels = voxelGrid.getOccupiedVoxels();
for (const auto& voxel : occupiedVoxels) {
    std::cout << "Voxel center: (" 
              << voxel.center.x << ", "
              << voxel.center.y << ", "
              << voxel.center.z << ")\n";
    std::cout << "Intersects " << voxel.atomIds.size() << " atom(s)\n";
}

// Access empty voxels (for visualization or other purposes)
const auto& emptyVoxels = voxelGrid.getEmptyVoxels();

// Get specific voxel by grid indices
const Voxel* voxel = voxelGrid.getVoxel(i, j, k);
if (voxel && voxel->occupied) {
    // Process occupied voxel
}
```

## Data Structures

### VoxelGrid

Contains the voxelized domain:
- `voxelSize_`: Edge length of each cubic voxel
- `dimensions_`: Number of voxels in x, y, z directions
- `occupiedVoxels_`: Vector of voxels that intersect with atoms
- `emptyVoxels_`: Vector of voxels that don't intersect with atoms
- `boundingBox_`: Bounding box used for voxelization

### Voxel

Represents a single voxel in the grid:
- `index`: Grid indices (i, j, k)
- `center`: Center point of the voxel
- `min`: Minimum corner coordinates
- `max`: Maximum corner coordinates
- `occupied`: Boolean flag for occupancy status
- `atomIds`: IDs of atoms that intersect this voxel

### HexMesh

Contains the generated mesh data:
- `nodes`: Vector of unique 3D node coordinates
- `elements`: Vector of element connectivity arrays (8 node indices per hexahedron)

### Node Ordering

Each hexahedral element uses the following node ordering:
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

- Nodes 0-3: Bottom face (z = min)
- Nodes 4-7: Top face (z = max)
- Node ordering follows right-hand rule

## Implementation Details

### Voxelization Algorithm

1. **Bounding Box Computation**: Calculate molecular bounding box with optional padding
2. **Grid Tessellation**: Divide bounding box into uniform cubic voxels based on specified edge length
3. **Occupancy Classification**: For each voxel, check if its center point lies within any atom sphere
4. **Voxel Storage**: Store occupied and empty voxels in separate vectors

### Occupancy Detection

The occupancy of a voxel is determined by checking if the voxel center intersects with any atom sphere:

```cpp
bool isOccupied = false;
for (each atom) {
    distance = ||voxel_center - atom_center||
    if (distance <= atom_radius) {
        isOccupied = true;
    }
}
```

This center-point method is:
- **Simple**: Single point test per voxel
- **Efficient**: O(n×m) where n = number of voxels, m = number of atoms
- **Conservative**: May miss small intersections at voxel corners
- **Adjustable**: Use smaller voxels for higher accuracy

### Mesh Generation

1. **Collect Occupied Voxels**: Process only voxels marked as occupied
2. **Compute Corner Nodes**: For each occupied voxel, compute its 8 corner coordinates
3. **Assign Unique Indices**: Use hash map to deduplicate shared nodes between adjacent voxels
4. **Build Connectivity**: Create element connectivity arrays mapping to unique node indices

The node deduplication uses a tolerance of 1e-12 for floating point coordinate comparison to handle numerical precision issues.

## Choosing Voxel Size

The voxel size parameter controls the mesh resolution:

| Voxel Size | Use Case | Characteristics |
|------------|----------|-----------------|
| 2.0 - 5.0 Å | Coarse mesh | Fast computation, low element count |
| 1.0 - 2.0 Å | Medium mesh | Balanced resolution and performance |
| 0.5 - 1.0 Å | Fine mesh | High resolution, captures small features |
| < 0.5 Å | Very fine mesh | Very high element count, slow for large molecules |

**Rule of thumb**: Voxel size should be comparable to or smaller than the smallest atomic radius (≈0.5 Å for hydrogen) for accurate representation.

## Performance Characteristics

- **Time Complexity**: O(n×m) where n = number of voxels, m = number of atoms
- **Space Complexity**: O(n) for storing voxels
- **Mesh Generation**: O(k) where k = number of occupied voxels

For a typical protein with 1000 atoms in a 30×30×30 Å box:
- 1 Å voxels: ~27,000 total voxels, typically 1-5% occupied
- 0.5 Å voxels: ~216,000 total voxels, typically 1-5% occupied

### OpenMP Parallelization

The mesh generation step is parallelized using OpenMP when available:

- **Parallelized Loop**: The computation of corner nodes for all occupied voxels runs in parallel
- **Thread Safety**: Uses indexed array assignment instead of push_back for thread-safe operations
- **Performance Impact**: Speedup is most noticeable for large structures with many occupied voxels (>10,000)
- **Configuration**: Control thread count with `OMP_NUM_THREADS` environment variable
- **Fallback**: Automatically runs serially if OpenMP is not available

Example:
```bash
# Use 8 threads for mesh generation
export OMP_NUM_THREADS=8
./voxel_demo large_protein.pdb 0.5
```

## Example Output

```
=== Voxel Grid Statistics ===
Voxel size: 1.00 Å
Grid dimensions: 10 x 13 x 11
Total voxels: 1430
Occupied voxels: 15 (1.0%)
Empty voxels: 1415 (99.0%)
Bounding box dimensions: 9.560 x 12.850 x 10.440 Å
Bounding box volume: 1282.51 Ų
=============================

Generated hexahedral mesh:
  Nodes: 91
  Elements: 15
  Node sharing efficiency: 24.2%
```

## Future Extensions

The separate storage of occupied and empty voxels enables:

- **Visualization**: Render occupied and empty regions differently
- **Multi-material Meshes**: Assign different material properties
- **Interface Tracking**: Identify boundary voxels between occupied and empty regions
- **Adaptive Refinement**: Further subdivide occupied voxels near boundaries
- **Image-Based Methods**: Integration with medical imaging data
- **Level Set Methods**: Distance field computation for implicit surfaces
