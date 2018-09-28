#include "timeev.hpp"
#include "setup.hpp"
#include <cmath>
#include <iostream>

template <class T>
double dd(const vector<T> &a, const vector<T> &b)
{
    complex<double> ret=0,c=0;
    for (unsigned i=0;i<a.size();i++)
    {
        c = a[i] - b[i];
        ret+=conj(c)*c;
    }
    return sqrt(ret.real());
}

double order(Setup &s, int dep, bool ubase, double dt)
{
    void (*cf)(Sparse_t&,vector<complex<double>>&,double,double,int);
    if (ubase){cf = &CFET;} else {cf = &CFET_nobase;}
    vector<complex<double>> v,w;
    v=s.v0();
    w=s.v0();
    cf(s.H, w, dt*0.5,0,dep);
    cf(s.H, w, dt*0.5,0,dep);
    cf(s.H, v, dt, 0, dep);
    double d = dd(v,w);
    return log10(d) + 10;
}

inline double arit_mean(double a, double b) { return (a+b)*0.5; }

double calibrate(Setup &s, int dep, bool ubase)
{
    double t_l = 1, t_h = 30, t, tt;
    do {
        t = arit_mean(t_l,t_h);
        tt = order(s,dep,ubase,t);
        if (tt < 0) { t_l = t; } else { t_h = t; }
        cout << t_h - t_l << endl;
    } while (t_h - t_l > 2e-1);
    return t_l;
}

void evolve(Setup &s, vector<complex<double>> &v, double t)
{
    static const int dep = 200;
    static Setup* init_ptr = nullptr;
    static void (*cf)(Sparse_t&,vector<complex<double>>&,double,double,int);
    static double d;

    Setup* now_ptr = &s;
    if (init_ptr != now_ptr) {
        init_ptr = now_ptr;
        cout << "Calibrating. " << flush;
        bool USE_BASE;
        USE_BASE = (s.b.vec.size() < 100000);
        if (USE_BASE){cf = &CFET;} else {cf = &CFET_nobase;}
        d = calibrate(s,dep,USE_BASE);
        cout << " Maximal timestep for Krylow depth " << dep << " is " << d << ". Saving base: " << USE_BASE << "\nDone.\n";
    }

    int sg = 2*(t>0.0) - 1;
    t = fabs(t);
    int n = floor( t/d );
    while (n--) { cf(s.H, v, d*sg, 0, dep); }
    cf(s.H, v, (t-floor(t/d)*d)*sg, 0, dep);
}
