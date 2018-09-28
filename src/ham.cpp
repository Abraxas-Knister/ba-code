#include "ham.hpp"

template <class T>
unsigned nofones(T x)
{
    unsigned n = 0;
    do n += x & 1; while ((x = x >> 1));
    return n;
}

template <class T>
inline bool oneat(const uint32_t &x, T i) { return x & (1 << i); }

Bookkeep::Bookkeep(uint8_t s, uint8_t p)
{
    cout << "Creating a Basis: " << flush;
    sites = s;
    particles = p;

    uint32_t j = ((1<<particles) - 1) << (sites - particles);
    const uint32_t last = (1 << particles) - 1;
    uint32_t i=0;
    do {
        if (nofones(j)==particles)
        {
            vec.push_back(j);
            ind.insert(pair<uint32_t,uint32_t>(j,i));
            i++;
        }
    } while (j-- != last);
    cout << "Done.\n";
}

int8_t CdC(uint8_t from, uint8_t to, uint32_t in, uint32_t &out, const Bookkeep &b) 
/* C_to^dagger after C_from. 
 *
 * :in: is number of basis element the CdC acts on
 * :out: is the basis element that is produced thereby
 * :b: is needed to calculate.
 *    - construct it by Bookkeep(sites_number, particles_number) beforehand
 * ----
 * returns:
 * :type: int
 * :val: 0, +/- 1
 *    - minus one power number of particles between from and to
 *    - zero if vector is destroyed by CdC whilst acting on it
 */
{
    uint32_t v = b.vec.at(in);
    if (!oneat(v,from)) return 0;
    if (from == to) return out=in, 1;
    if (oneat(v,to)) return 0;
    out = b.ind.at((v | (1<<to)) & (~(1<<from)));

    int c = from < to;

    for (unsigned i=0; i<b.sites; i++)
    {
        if (oneat(v,i)) c += (i<from) + (i<to);
    }
    if (c&1) return -1;
    return 1;  
}

void genham(Sparse_t &hammy, const vector<vector<int8_t>> &matham, const Bookkeep &b, const bool ferm)
{
    cout << "Creating a Hamiltonian: " << flush;
    hammy.clear();
    
    unsigned dim = b.vec.size();
    uint64_t ptr = 0;

    cout << "The Dimension is " << dim << ".\n";
    int c=0; // c is just for the .... loading bar
    hammy.row_ptr.push_back(ptr); 
    for (uint32_t in = 0; in < dim; in++)
    {
        #pragma omp parallel for
        for (uint8_t from = 0; from < b.sites; from++)
        {
            for (uint8_t to = 0; to < b.sites; to++)
            {
                int8_t h = matham.at(from).at(to);
                if (h)
                {
                    uint32_t out;
                    int8_t sg = CdC(from,to,in,out,b);
                    if (sg)
                    {
                        #pragma omp critical
                        {
                            if (!ferm) {sg = 1;}
                            hammy.val.push_back(sg*h);
                            hammy.col_ind.push_back(out);
                            ptr++;
                        }
                    } else {
                        continue;
                    }
                }
            }
        }
        hammy.row_ptr.push_back(ptr); 
        // Loading bar ........................................
        if (in * 100 / dim > c) {c++;cout << "." << flush;}
    }
    cout << "\nDone.\n";
}
