#include "sparse.hpp"
#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;
using namespace BA;

Sparse::Sparse()
{
    col_ind = {};
    val = {};
    row_ptr = {};
}

void Sparse::xx(const vector<complex<double>> &x, vector<complex<double>> &w) const
{
    uint32_t dim = row_ptr.size() - 1;
    w.assign(dim,0);
    #pragma omp parallel for
    for (uint32_t i=0; i<dim; i++)
    {
        for (uint32_t j= row_ptr[i]; j<row_ptr[i+1]; j++)
        {
            w[i] += ((complex<double>) val[j])*x[col_ind[j]];
        }
    }
}

void Sparse::clear()
{
    col_ind.clear();
    val.clear();
    row_ptr.clear();
}
