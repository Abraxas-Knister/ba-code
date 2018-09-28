#include "timeev.hpp"

void toMat(MatrixXcd&, const vector<double>&, const vector<double>&);
void expdotv(Sparse_t&, vector<complex<double>>&, double, int);

void CFET(Sparse_t &H, vector<complex<double>> &v, double dt, double t, int depth)
{
    expdotv(H,v,dt,depth); // H assumed const
///  - - - - - -  CFET needs this if H = H(t)
/*
 *    double x[] = {
 *      -0.387298334620741688517926539978,
 *       0.5,
 *       0.387298334620741688517926539978
 *   };
 *
 *   double g[] = {
 *       0.005776500145309697885851900391,
 *      -0.033333333333333333333333333333,
 *       0.302556833188023635447481432942,
 *      -0.027777777777777777777777777778,
 *       0.511111111111111111111111111111
 *   };
 *  
 */

////////////
}


void toMat(MatrixXcd &mat, const vector<double> &a, const vector<double> &b)
{
    int dim = a.size();
    for (int i=0;i<dim;i++)
    {
        mat(i,i) = a[i];
        if (i) {
            mat(i,i-1) = b[i];
            mat(i-1,i) = b[i];        
        }
    }
}

void expdotv(Sparse_t &H, vector<complex<double>> &v, double mult, int depth)
{
    vector<double> a,b;
    vector<vector<complex<double>>> base;
    lanczos(H,v,a,b,base,depth);

    int dim = a.size();
    MatrixXcd H_mat = MatrixXcd::Zero(dim,dim);
    toMat(H_mat,a,b);

    H_mat *= mult/1i;
    VectorXcd u = VectorXcd::Zero(dim);
    VectorXcd tmp = VectorXcd::Zero(dim);
    tmp(0)=1;
    u = H_mat.exp()*tmp;

    int bigdim = v.size();
    v.clear();
    v.resize(bigdim,0);

    for (int j=0;j<bigdim;j++)
    {
        for (int i=0;i<dim;i++)
        {
            v[j] += u(i) * base[i][j];
        }
    }
}

void CFET_nobase(Sparse_t &H, vector<complex<double>> &v, double dt, double t, int depth)
{
    vector<complex<double>> v_bak=v;
    vector<double> a,b;
    lanczos_coeff(H,v,a,b,depth);
    MatrixXcd H_mat = MatrixXcd::Zero(depth,depth);
    toMat(H_mat,a,b);

    H_mat *= dt/1i;
    VectorXcd u = VectorXcd::Zero(depth);
    VectorXcd tmp = u;
    tmp(0)=1;
    u = H_mat.exp()*tmp;

    v = v_bak;
    lanczos_assemble(H,v,a,b,u);
}
