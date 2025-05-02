#ifndef DDM_H
#define DDM_H

#include <vector>
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>

#include "io_handler.h"
#include "matrix_tools.h"
#include "global_defs.h"
#include "tools.h"

namespace GroupIP
{
    std::vector<std::vector<int>> vertex_hspaces_adjacency(const Matrix &A, const Vector &b);
    dd_MatrixPtr generate_matrix_for_triangulation(const std::vector<Vector> &cone_rays, int dim);
    void triangulation(const Matrix &A, std::vector<std::vector<int>> &simple_cones, const std::vector<Vector> &cone_rays, int dim);
    std::vector<std::vector<int>> get_cones(Matrix &A, Vector &b);
}

#endif // DDM_H