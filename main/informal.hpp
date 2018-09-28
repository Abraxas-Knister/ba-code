#include "timeev.hpp"
#include "setup.hpp"
#include "density.hpp"
/*
 * This header does things that are not discussed in the thesis 
 *   but can come in handy if one would like to discuss fermions
 *   in a rather numerous number of degrees of freedom
 * we are solving heisenberg eom for the density operator here.
 */

#include <cmath>
#include <iomanip>

int main(int argc,char* argv[])
{
    // argc isn't used -> there will be an error that ist displaying it
    string tmp = argv[1];
    Setup s(tmp);
    // get access to the setup you configured using GUI.py
    /*
     * this setup information is
     * :b: struct Bookkeep
     *   necceccary information about the hilbert space
     * :H: Sparse CRS
     *   built up from the adjacency matrix
     * :ferm: bool
     *   spinless fermions if true otherwise hard-core bosons
     * :subs: vector
     *   indices of subsystem places (red in GUI.py)
     * :v0: function
     *   call to get initial config
     * :noins: function
     *   call with argument vector to get 'no in subsystem'
     */
    cout << "Done building the setup. Accumulating data for " + "information purposes\n";
    ofstream f(tmp+"info");
    // setup data ofstream
    vector<complex<double>> v=s.v0(), w;
    // initial configuration, temporary storage
    double dt = 1;
    for (unsigned k=1; k<=50; k++)
    {
        evolve(s,v,dt); w=v;
        /*
         * five routines for time evolution are there:
         *  CFET / CFET_nobase :
         *   it was planned to set up a time evolution that can handle a time dependant Hamiltonian
         *     compute -> exp(- \int H dt)
         *     CFETs are things that can do it nicely, but here its just a function name
         *   difference: one stores the  lanczos base one doesnt
         *  
         *  Ev_future / Ev_past :
         *   Call CFET n times either to go into the past or to go into the future
         *  evolve :
         *   intelligently choose a bigest time intervall 
         *     that can be evolved over without accuracy problems
         *   compare your timestep with the computed one
         *     go as many with the big one and the rest with a small one
         *   (this thing can also go into the past)  
         */
        disturbe_orth(w, s.b, 0.3);
        /*
         * three routines for disturbement:
         *  disturbe: just apply phasekick
         *  disturbe_orth: also orthonormalize
         *  disturbe_op: do it in limit of small time
         * first two are dependant of some parameter alpha
         */
        evolve(s,w,-k*dt);
        f << setprecision(15) << dt*k;
        f << " " << setprecision(15) << s.noins(v)/s.subs.size(); // <O>_v at time t in units of subsystem size
        f << " " << setprecision(15) << densN(t,"20")/s.subs.size(); // <0>(t)  if it were fermions
        // densN comes from the density.hpp
        f << " " << setprecision(15) << s.noins(w) // <O> using ket U(-t)PU(t)v
        f << endl;
        f.flush();
    }
    f.close();
    cout << "Done."
}
