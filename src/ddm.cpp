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
                    hspaces_adjacency[i].push_back(elem - 1);
            }
        }

        dd_FreeMatrix(dd_A);
        dd_FreeMatrix(G);
        dd_FreePolyhedra(poly);
        dd_free_global_constants();

        return hspaces_adjacency;
    }

    void polytope_preprocessing(Matrix &A, Vector &b)
    {
        dd_ErrorType err;
        dd_MatrixPtr dd_A = get_cdd_system(A, b);

        dd_rowset impl_lin, redset;
        dd_rowindex newpos;
        dd_MatrixCanonicalize(&dd_A, &impl_lin, &redset, &newpos, &err);

        std::vector<int> linearity_row_idxs;
        std::vector<int> inequality_row_idxs;
        for (int i = 0; i < dd_A->rowsize; ++i)
        {
            if (set_member(i + 1, dd_A->linset))
            {
                linearity_row_idxs.push_back(i);
                continue;
            }
            inequality_row_idxs.push_back(i);
        }

        if (linearity_row_idxs.size() == 0)
        {
            return;
        }

        Matrix E;
        Vector Eb;
        for (auto i : linearity_row_idxs)
        {
            Eb.push_back(*dd_A->matrix[i][0]);
            Vector row;
            for (int j = 1; j < dd_A->colsize; ++j)
            {
                row.push_back(-(*dd_A->matrix[i][j]));
            }
            E.push_back(row);
        }

        Matrix new_A;
        Vector new_b;

        for (auto i : inequality_row_idxs)
        {
            Vector row;
            new_b.push_back(*dd_A->matrix[i][0]);
            for (int j = 1; j < dd_A->colsize; ++j)
            {
                row.push_back(-(*dd_A->matrix[i][j]));
            }
            new_A.push_back(row);
        }

        Matrix H, U;
        hermite_normal_form(E, H, U);

        // solution in new variables
        Vector new_vars_solution(H[0].size());
        if (Eb[0] % H[0][0] != 0)
        {
            std::cout << "Num integer points: 0\n";
            exit(1);
        }

        new_vars_solution[0] = Eb[0] / H[0][0];

        for (int i = 1; i < H.size(); ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                Eb[i] -= H[i][j] * new_vars_solution[j];
            }

            if (Eb[i] % H[i][i] != 0)
            {
                std::cout << "Num integer points: 0\n";
                exit(1);
            }

            new_vars_solution[i] = Eb[i] / H[i][i];
        }

        Vector constant_terms(A[0].size(), 0);
        int num_vars = E[0].size();
        int num_new_vars = E[0].size() - E.size();
        int start_new_var_idx = E.size();
        Matrix new_AA(inequality_row_idxs.size(), Vector(num_new_vars, 0));
        Vector new_bb(inequality_row_idxs.size(), 0);

        for (int i = 0; i < num_vars; ++i)
        {
            for (int j = 0; j < E.size(); ++j)
            {
                constant_terms[i] += U[i][j] * new_vars_solution[j];
            }
        }

        for (int m = 0; m < new_AA.size(); ++m)
        {
            for (int i = 0; i < num_new_vars; ++i)
            {
                for (int j = 0; j < num_vars; ++j)
                {
                    new_AA[m][i] += new_A[m][j] * U[j][i + start_new_var_idx];
                }
            }
        }

        for (int m = 0; m < new_A.size(); ++m)
        {
            new_bb[m] = new_b[m];
            for (int i = 0; i < num_vars; ++i)
            {
                new_bb[m] -= new_A[m][i] * constant_terms[i];
            }
        }

        A = new_AA;
        b = new_bb;
    }

    std::vector<std::vector<int>> get_cones(Matrix &A, Vector &b)
    {
        polytope_preprocessing(A, b);
        int dim = A[0].size();

        std::vector<std::vector<int>> simple_cones;
        auto cones = vertex_hspaces_adjacency(A, b);

        for (auto cone : cones)
        {
            if (cone.size() == dim)
            {
                simple_cones.push_back(cone);
            }
            else
            {
                Matrix Acone;
                for (auto normal_idx : cone)
                {
                    Acone.push_back(A[normal_idx]);
                }
                triangulation(A, simple_cones, Acone, dim);
            }
        }

        return simple_cones;
    }

    dd_MatrixPtr generate_matrix_for_triangulation(const std::vector<Vector> &cone_rays, int dim)
    {
        dd_set_global_constants();
        int num_rays = cone_rays.size();
        dd_MatrixPtr matrix = dd_CreateMatrix(num_rays + 1, dim + 2);
        matrix->numbtype = dd_Rational;
        matrix->representation = dd_Generator;

        mpq_class x;
        for (int i = 0; i < num_rays; ++i)
        {
            for (int j = 0; j < dim; j++)
            {
                dd_set_si(matrix->matrix[i][j + 2], cone_rays[i][j].get_si());
            }
        }
        return matrix;
    }

    void triangulation(const Matrix &A, std::vector<std::vector<int>> &simple_cones, const std::vector<Vector> &cone_rays, int dim)
    {
        dd_MatrixPtr matrix = generate_matrix_for_triangulation(cone_rays, dim);
        int num_rays = cone_rays.size();

        dd_set_si(matrix->matrix[num_rays][1], 1);

        // compute heights
        for (int i = 0; i < num_rays; ++i)
        {
            auto height = uniform_random_number(1, 100);
            dd_set_si(matrix->matrix[i][1], height);
            dd_set_si(matrix->matrix[i][0], 1);
        }

        // compute triangulation
        dd_ErrorType error;
        dd_PolyhedraPtr poly = dd_DDMatrix2Poly(matrix, &error);
        dd_MatrixPtr inequalities = dd_CopyInequalities(poly);
        dd_SetFamilyPtr incidence = dd_CopyIncidence(poly);

        int num_inequalities = inequalities->rowsize;
        for (int i = 0; i < num_inequalities; i++)
        {
            if (!set_member(i + 1, inequalities->linset) && !set_member(num_rays + 1, incidence->set[i]))
            {
                int n_rays = set_card(incidence->set[i]);
                auto ray_set = incidence->set[i];
                if (n_rays == dim)
                {
                    std::vector<int> simple_cone;
                    for (int j = 1; j < num_rays + 1; ++j)
                    {
                        if (set_member(j, ray_set))
                        {
                            int idx = std::find(A.begin(), A.end(), cone_rays[j - 1]) - A.begin();
                            simple_cone.push_back(idx);
                        }
                    }
                    simple_cones.push_back(simple_cone);
                }
                else
                {
                    Matrix Acone;
                    for (int j = 1; j < num_rays + 1; ++j)
                    {
                        if (set_member(j, ray_set))
                        {
                            Acone.push_back(cone_rays[j - 1]);
                        }
                    }
                    triangulation(A, simple_cones, Acone, dim);
                }
            }
        }
    }
}
