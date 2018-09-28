#pragma once
#include "ham.hpp"
#include <iostream>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
typedef BA::Sparse Sparse_t;

struct Setup {
    Bookkeep b;
    Sparse_t H;
    bool ferm;
    vector<uint8_t> subs;
    vector<complex<double>> v0(void);
    Setup(string);
    double noins(const vector<complex<double>>&);
};
    
void load(string, vector<vector<int8_t>>&, vector<uint8_t>&, uint8_t&, bool&);
