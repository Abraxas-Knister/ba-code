#pragma once
#include <iostream>
#include <cstdint>
#include <vector>
#include <omp.h>
#include <map>
#include "sparse.hpp"

using namespace std;
typedef BA::Sparse Sparse_t;

struct Bookkeep {
    vector<uint32_t> vec;
    map<uint32_t,uint32_t> ind;
    
    uint8_t sites;
    uint8_t particles;
    Bookkeep(){}
    Bookkeep(uint8_t ,uint8_t);
};

void genham(Sparse_t&, const vector<vector<int8_t>>&, const Bookkeep&, const bool);
