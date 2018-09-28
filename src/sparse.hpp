#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cstdint>

using namespace std;

namespace BA {
struct Sparse
{
    vector<uint32_t> col_ind;
    vector<uint64_t> row_ptr;
    vector<int8_t> val;
    Sparse();
    void xx(const vector<complex<double>> &, vector<complex<double>> &) const;
    void clear();
    void print();
};}
