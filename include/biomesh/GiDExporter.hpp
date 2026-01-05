#pragma once

#include "biomesh/HexMesh.hpp"
#include <string>

namespace biomesh {

/**
 * @brief GiD mesh file exporter
 * 
 * Exports hexahedral meshes to GiD proprietary format (.msh files).
 * GiD is a universal pre/post-processing software widely used in FEM/CFD.
 * 
 * The format uses 1-based indexing and follows the structure:
 * - MESH header with dimension and element type
 * - Coordinates section with node positions
 * - Elements section with connectivity
 */
class GiDExporter {
public:
    /**
     * @brief Export hexahedral mesh to GiD format (.msh file)
     * @param mesh The hexahedral mesh to export
     * @param filename Output file path (should end with .msh)
     * @return true if export was successful, false otherwise
     */
    static bool exportToGiD(const HexMesh& mesh, const std::string& filename);
};

} // namespace biomesh
