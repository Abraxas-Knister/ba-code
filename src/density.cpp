#include "density.hpp"
#include "setup.hpp"
#include <cstdint>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

double trS(const MatrixXcd &p, const vector<uint8_t> a = {})
{
    static bool FF = true;
    static vector<uint8_t> subs;
    if (FF)
    {
        FF = false;
        subs = a;
        return 0;
    }
    complex<double> ret=0;
    for (auto i:subs) ret += p(i,i);
    return fabs(ret);
}

void init(string proj, MatrixXcd &H, MatrixXcd &p)
{
    vector<vector<int8_t>> H_m;
    vector<uint8_t> subs;
    bool f;
    uint8_t ps;
    load(proj,H_m,subs,ps,f);
    H.resize(H_m.size(),H_m.size());
    for (unsigned i=0;i<H_m.size();i++)
    { 
        for (unsigned j=0;j<H_m.size();j++)
        {
            H(i,j) = H_m[i][j];
        }
    }

    p = MatrixXcd::Zero(H_m.size(),H_m.size());
    for (auto i:subs) {p(i,i) = 1;}
    trS(p,subs);
}


double densN(double t, string proj)
{
    static MatrixXcd H,p;
    static string pr = "";
    if (pr != proj)
    {
        pr=proj;
        init(proj,H,p);
    }
    MatrixXcd tmp = ((t/1i*H).exp())*p*((t*1i*H).exp());
    return trS(tmp);
}       


