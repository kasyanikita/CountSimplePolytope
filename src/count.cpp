#include <iostream>
#include <fstream>

#include "counter.h"
#include "io_handler.h"
#include "ddm.h"

using namespace GroupIP;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " filepath\n";
        std::cerr << "Error: Expected 2 arguments, but got " << argc - 1 << ".\n";
        return 1;
    }

    std::string filepath = argv[1];
    Matrix A;
    Vector b;

    read_data(A, b, filepath);

    int dim = A[0].size();

    auto hspace_adjacency = vertex_hspaces_adjacency(A, b);

    for (auto cone_ids : hspace_adjacency)
    {
        if (cone_ids.size() != dim)
        {
            std::cerr << "The polytope isn't simple!\n";
            std::exit(1);
        }
    }

    std::vector<Matrix> A_cones;
    std::vector<Vector> b_cones;

    for (auto cone_idxs : hspace_adjacency)
    {
        Matrix A_cone;
        Vector b_cone;
        for (auto idx : cone_idxs)
        {
            A_cone.push_back(A[idx - 1]);
            b_cone.push_back(b[idx - 1]);
        }
        A_cones.push_back(A_cone);
        b_cones.push_back(b_cone);
    }

    HyperplaneAvoidSolver hyperplane_avoid_vector(A);
    auto c = hyperplane_avoid_vector.get_vector(A_cones);

    mpq_class res = 0;

    for (int_t i = 0; i < A_cones.size(); ++i)
    {
        res += cone_evaluation(A_cones[i], b_cones[i], c);
    }

    std::cout << "Number of integer points: " << res << std::endl;

    return 0;
}