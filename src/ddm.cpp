#include "ddm.h"

namespace GroupIP
{
    std::vector<std::vector<int>> vertex_hspaces_adjacency(const Matrix &A, const Vector &b)
    {
        dd_PolyhedraPtr poly;
        dd_MatrixPtr dd_A, G;
        dd_SetFamilyPtr GI;
        dd_rowrange m;
        dd_colrange d;
        dd_ErrorType err;

        dd_set_global_constants();

        assert(A.size() > 0);

        m = A.size();
        d = A[0].size() + 1;
        dd_A = dd_CreateMatrix(m, d);

        for (int i = 0; i < m; ++i)
        {
            dd_set_si(dd_A->matrix[i][0], b[i].get_si());
            for (int j = 1; j < d; ++j)
            {
                dd_set_si(dd_A->matrix[i][j], -A[i][j - 1].get_si());
            }
        }

        dd_A->representation = dd_Inequality;
        poly = dd_DDMatrix2Poly(dd_A, &err);

        G = dd_CopyGenerators(poly);
        GI = dd_CopyIncidence(poly);

        int rowsize = G->rowsize;
        int colsize = G->colsize;

        if (rowsize == 0)
        {
            std::cerr << "Empty polyhedron!\n";
            std::exit(1);
        }

        std::vector<std::vector<double>> vertices(
            rowsize, std::vector<double>(colsize - 1, 0));

        for (int i = 0; i < G->rowsize; ++i)
        {
            if (*G->matrix[i][0] == 0)
            {
                std::cerr << "The polyhedron is unbounded!\n";
                std::exit(1);
            }
            for (int j = 0; j < G->colsize - 1; ++j)
            {
                vertices[i][j] = *G->matrix[i][j + 1];
            }
        }

        std::vector<std::vector<int>> hspaces_adjacency(vertices.size(),
                                                        std::vector<int>());
        for (int i = 0; i < GI->famsize; ++i)
        {
            long elem;

            for (elem = 1; elem <= GI->set[i][0]; elem++)
            {
                if (set_member(elem, GI->set[i]))
                    hspaces_adjacency[i].push_back(elem);
            }
        }

        dd_FreeMatrix(dd_A);
        dd_FreeMatrix(G);
        dd_FreePolyhedra(poly);
        dd_free_global_constants();

        return hspaces_adjacency;
    }
}