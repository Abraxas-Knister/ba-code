#pragma once
#define SEED 928734
#include "sparse.hpp"
#include "ham.hpp"
#include <complex>
#include <cstdlib>
#include <unsupported/Eigen/MatrixFunctions>

struct Setup;

using namespace std;
using namespace Eigen;
using namespace complex_literals;
typedef BA::Sparse Sparse_t;

complex<double> dot(const vector<complex<double>>&, const vector<complex<double>>&);

void CFET(Sparse_t&, vector<complex<double>>&, double,double,int); 
void lanczos(const Sparse_t&, vector<complex<double>>&, vector<double>&, vector<double>&, vector<vector<complex<double>>>&, int);

void CFET_nobase(Sparse_t&, vector<complex<double>>&, double,double,int); 
void lanczos_coeff(const Sparse_t&, vector<complex<double>>&, vector<double>&, vector<double>&, int);
void lanczos_assemble(const Sparse_t&, vector<complex<double>>&, const vector<double>&, const vector<double>&, const VectorXcd&);

void evolve(Setup &, vector<complex<double>> &, double);
void Ev_future(Sparse_t&, vector<complex<double>>&, double,int,int);
void Ev_past(Sparse_t&, vector<complex<double>>&, double,int,int);
void disturbe(vector<complex<double>>&,const Bookkeep&,double,int=SEED);
void disturbe_orth(vector<complex<double>>&,const Bookkeep&,double,int=SEED);
void disturbe_op(vector<complex<double>>&,const Bookkeep&,int=SEED);
