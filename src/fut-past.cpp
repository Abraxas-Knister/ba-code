#include "timeev.hpp"

void Ev_future(Sparse_t &H, vector<complex<double>> &v, double dt, int steps, int depth)
{
    for (int i=0;i<steps;i++)
    {
        CFET(H,v,dt,dt*i, depth);
    }
}

void Ev_past(Sparse_t &H, vector<complex<double>> &v, double dt, int steps, int depth)
{
    for (int i=steps;i>0;i--)
    {
        CFET(H,v,-dt,dt*i, depth);
    }
}
