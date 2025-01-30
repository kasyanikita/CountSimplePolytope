#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#include <iostream>
#include <fstream>
#include <string>

#include "global_defs.h"

namespace GroupIP
{
    void read_data(Matrix &A, Vector &b, const std::string &filepath);
    void print_matrix(const Matrix &matrix);
    void print_vector(const Vector &vector);
}

#endif // IO_HANDLER_H