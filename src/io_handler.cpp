#include "io_handler.h"

namespace GroupIP
{
    void read_data(Matrix &A, Vector &b, const std::string &filepath)
    {
        std::fstream file(filepath);

        if (!file.is_open())
        {
            std::cerr << "Failed to open the file " << filepath << std::endl;
            std::exit(1);
        }

        int n_rows = 0;
        int n_cols = 0;
        file >> n_rows >> n_cols;

        for (int i = 0; i < n_rows; ++i)
        {
            Vector row;
            for (int j = 0; j < n_cols; ++j)
            {
                int value;
                file >> value;
                if (j == 0)
                {
                    b.push_back(value);
                }
                else
                {
                    row.push_back(-value);
                }
            }
            A.push_back(row);
        }
    }

    void print_matrix(const Matrix &matrix)
    {
        int n_rows = matrix.size();

        for (int i = 0; i < n_rows; ++i)
        {
            print_vector(matrix[i]);
            std::cout << std::endl;
        }
    }

    void print_vector(const Vector &vector)
    {
        int vec_size = vector.size();

        for (int i = 0; i < vec_size; ++i)
        {
            std::cout << vector[i] << " ";
        }
    }
}