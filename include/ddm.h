#ifndef DDM_H
#define DDM_H

#include <vector>
#include <cddlib/setoper.h>
#include <cddlib/cdd.h>

#include "global_defs.h"

namespace GroupIP
{
    std::vector<std::vector<int>> vertex_hspaces_adjacency(const Matrix &A, const Vector &b);
}

#endif // DDM_H