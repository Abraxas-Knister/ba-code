#include "timeev.hpp"

complex<double> dot(const vector<complex<double>> &a, const vector<complex<double>> &b)
{
    complex<double> ret;
    for (unsigned i = 0; i<b.size(); i++) ret += conj(a.at(i))*b.at(i);
    return ret;
}

void scale(complex<double> mult, vector<complex<double>> &v)
{
    for (unsigned i=0;i<v.size();i++) v[i] *= mult;
}

void axpy(complex<double> m, const vector<complex<double>> &v, vector<complex<double>> &w)
{
    for (unsigned i = 0; i<w.size(); i++)
    {
        w[i] += m*v[i];
    }
}

double norm(const vector<complex<double>> &v) { return sqrt(fabs(dot(v,v))); }

void lanczos(const Sparse_t &H, vector<complex<double>> &v, vector<double> &a, vector<double> &b, vector<vector<complex<double>>> &base, int depth)
{
    a.clear(), b.clear(), base.clear();
    b.push_back(norm(v));
    scale(1/b[0], v);
    base.push_back(v);

    vector<complex<double>> w,tmp;
    H.xx(v,w);
    a.push_back(dot(v,w).real());
    axpy(-a[0],v,w);
    b.push_back(norm(w));

    int n = 1;
    while (fabs(b[n]) > 1e-10)
    {
        scale(1/b[n], w);
        base.push_back(w);

        scale(-b[n], v);
        v.swap(w);
        
        H.xx(v,tmp);
        for (unsigned i=0;i<w.size();i++) w[i] += tmp[i];

        a.push_back(dot(v,w).real());
        axpy(-a[n],v,w);
        b.push_back(norm(w));
        n++;
        if (n>=depth) return;
    }
}

void lanczos_coeff(const Sparse_t &H, vector<complex<double>> &v, vector<double> &a, vector<double> &b, int depth)
{
    a.clear(), b.clear();
    b.push_back(norm(v));
    scale(1/b[0], v);
    vector<complex<double>> w,tmp;
    H.xx(v,w);
    a.push_back(dot(v,w).real());
    axpy(-a[0],v,w);
    b.push_back(norm(w));

    int n = 1;
    while (fabs(b[n]) > 1e-10)
    {
        scale(1/b[n], w);
        scale(-b[n], v);
        v.swap(w);
        H.xx(v,tmp);
        for (unsigned i=0;i<w.size();i++) w[i] += tmp[i];
        a.push_back(dot(v,w).real());
        axpy(-a[n],v,w);
        b.push_back(norm(w));
        n++;
        if (n>=depth) return;
    }
}

void lanczos_assemble(const Sparse_t &H, vector<complex<double>> &v, const vector<double> &a, const vector<double> &b, const VectorXcd &u)
{
    vector<complex<double>> w,tmp,ret(v.size(),0);
    scale(1/b[0], v);
    for (unsigned l=0;l<ret.size();l++) { ret[l] +=u(0)*v[l]; }
    H.xx(v,w);
    axpy(-a[0],v,w);
    int n = 1;
    for (unsigned k=1; k<u.size(); k++)
    {
        scale(1/b[n], w);
        for (unsigned l=0;l<ret.size();l++) { ret[l] +=u(k)*w[l]; }
        scale(-b[n], v);
        v.swap(w);
        H.xx(v,tmp);
        for (unsigned i=0;i<w.size();i++) w[i] += tmp[i];
        axpy(-a[n],v,w);
        n++;
    }
    v.swap(ret);
}
