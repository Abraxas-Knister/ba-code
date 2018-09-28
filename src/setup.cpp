#include "setup.hpp"
#include "timeev.hpp" // bc. of use of dot(v,w) in noins.

typedef istringstream ins;
typedef ostringstream outs;

uint8_t num(uint16_t);

Setup::Setup(string proj)
{
    vector<vector<int8_t>> tmp_ham;
    uint8_t p;
    load(proj,tmp_ham,subs,p,ferm);
    b = Bookkeep(tmp_ham.size(),p);
    genham(H,tmp_ham,b,ferm);
    noins(v0());
}

double Setup::noins(const vector<complex<double>> &v)
{
    static vector<double> mult;
    static bool FF = true;
    if (FF)
    {
        FF = false;
        cout << "Assigning multiplicators for the subs particle number op: " << flush;
        mult.assign(b.vec.size(),0);
        #pragma omp parallel for
        for (unsigned k=0;k<b.vec.size();k++)
        {
            int c=0;
            for (auto i:subs)
            {
                if (b.vec[k] & (1<<i)) c++;
            }
            mult[k] = c;
        }
        cout << "Done.\n";
        return 0;
    }

    complex<double> ret=0;
    for (unsigned k=0; k<b.vec.size(); k++)
    {
        ret += conj(v[k])*v[k]*mult[k];
    }

    ret /= dot(v,v);
    return ret.real();
}

void load(string project, vector<vector<int8_t>> &H, vector<uint8_t> &subs, uint8_t &particles, bool &ferm)
{
    string s;
    int i;
    ifstream valsrc(".projects/" + project + "_vars");
    getline(valsrc, s);
    (ins) s >> i;
    particles = i;
    getline(valsrc, s);
    (ins) s >> i;
    ferm = i;
    subs.clear();
    while (valsrc.good())
    {
        getline(valsrc, s);
        (ins) s >> i;
        subs.push_back(i);
    }

    ifstream matsrc(".projects/" + project + "_nt");
    vector<int8_t> coeff;
    while (matsrc.good())
    {
        getline(matsrc, s);
        (ins) s >> i;
        coeff.push_back(i);
    }

    uint8_t dim = num(coeff.size());
    H.clear();
    for (uint8_t i=0;i<dim;i++) H.push_back(vector<int8_t>(dim,0));
    uint8_t pref = 0;
    for (uint8_t i=dim; i>0; i--)
    {
        for (uint8_t j=0; j<i; j++)
        {
            H[j][j+dim-i] = coeff[pref + j];
            H[j+dim-i][j] = coeff[pref + j];
        }
        pref += i;
    }
}
     
uint8_t num(uint16_t n)
/*
 *   :n: -- some positive number
 * :ret: -- sum(i=1,ret)i shall equal n
 *           if this is not possible -1 is returned
 */
{
    uint8_t ret = 0;
    while (1)
    {
        if (n==0) return ret;
        if (n<0) return -1;
        n -= ++ret;
    }
}

vector<complex<double>> Setup::v0()
{
    if (subs.size() != b.particles) {cout << "Don't come up with this Shit!\n We're doing only fully filled subsystems right now!\n"; throw;}
    else {
        vector<complex<double>> v;
        uint32_t w=0;
        for (auto i: subs) w |= 1<<i;
        uint32_t i = b.ind.at(w);
        v.assign(b.vec.size(),0);
        v[i]=1;
        return v;
    }
}

