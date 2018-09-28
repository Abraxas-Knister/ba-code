#include "timeev.hpp"
#include <cmath>

void _init(vector<double> &ww, const Bookkeep &b, int seed = SEED)
{
    static vector<double> w(b.vec.size(),0);
    static int prev=0;
    if (prev!=seed)
    {
        prev=seed;
        srand(seed);
        cout << "Newly creating random numbers\n";
        vector<double> r(b.sites, 0); 
        for (unsigned k=0; k<r.size(); k++)
        {
            r[k] = ((double)rand())/
                (1.0 + (double)RAND_MAX);
        }
        for (unsigned k=0; k<b.vec.size(); k++)
        {
            for (unsigned l=0;l<b.sites;l++)
            {
                if (b.vec[k] & 1<<l) w[k] += r[l];
            }
        }
        for (unsigned k=0;k<w.size();k++) w[k] /= sqrt((double) b.sites);
    }
    ww = w;
}

void disturbe(vector<complex<double>> &v, const Bookkeep &b, double strength,int s)
{
    vector<double> w;
    _init(w,b,s);
    for (unsigned k=0;k<v.size();k++)
        v[k] *= exp(-1i*w[k]*strength); 
}

void disturbe_orth(vector<complex<double>> &v, const Bookkeep &b, double strength,int s)
{
    vector<double> w;
    _init(w,b,s);
    vector<complex<double>> zz = v;
    for (unsigned k=0;k<v.size();k++)
        zz[k] *= exp(-1i*w[k]*strength); 
    complex<double> c = dot(v,zz);
    for (unsigned k=0;k<v.size();k++)
        v[k] = zz[k] - v[k]*c;
    c = sqrt(dot(v,v).real());
    for (unsigned k=0;k<v.size();k++)
        v[k]/=c;
}

void disturbe_op(vector<complex<double>> &v, const Bookkeep &b, int s)
{
    vector<double> w;
    _init(w,b,s);
    vector<complex<double>> d=v;
    for (unsigned k=0;k<v.size();k++)
        d[k] *=w[k];
    complex<double> c=dot(v,d);
    for (unsigned k=0;k<v.size();k++)
        v[k] = d[k] - c*v[k];
    c = sqrt(dot(v,v).real());
    for (unsigned k=0;k<v.size();k++)
        v[k]/=c;
}

