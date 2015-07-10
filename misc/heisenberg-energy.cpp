// -*- C++ -*- $Id$

// Bethe-ansatz energies for the Heisenberg model for finite PBC

// From arxiv:cond-mat/9809163

// The solution is specified by a set of numbers {z_r}, which 
// are solutions of the equation
// N \phi(z_i) = 2 pi I_i + \sum_{j \neq i} [(z_i - z_j)/2]
// for i = 0, .... , r-1
// and I_i are the 'bethe ansatz quantum numbers'.  For the groundstate,
// these are simply I_i = -N/4 + 0.5*(i+1)
// and phiz(z) = 2 arctan z is the 'bare momentum' of the magnon.
// The energy of each magnon ie epsilon(z) = dk/dz = -2/(1 + z^2)
// Note that generally the {z_r} are complex, but in this case
// they are real.

#include "common/trace.h"
#include "common/math_const.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>

using namespace std;

// Bethe quantum numbers
double Ir(int N, int r)
{
   DEBUG_CHECK(r >= 0 && r < N);
   return -0.25*N + 0.5 + r;
}

double phi(double z)
{
   return 2 * atan(z);
}

// iterative scheme

void Iterate(int N, std::vector<double>& New_z, std::vector<double> const& z)
{
   int const rMax = z.size();
   for (int i = 0; i < rMax; ++i)
   {
      double Sum = 0.0;
      for (int j = 0; j < rMax; ++j)
      {
	 if (j != i)
	 {
	    Sum += phi((z[i] - z[j])/2);
	 }
      }
      New_z[i] = tan((math_const::pi * Ir(N,i) + 0.5*Sum) / N);
   }
}

// magnon energy
double epsilon(double z)
{
   return -2.0 / (1 + z*z);
}

int main(int argc, char** argv)
{
   if (argc < 2)
   {
      cerr << "Calculate the exact energy for the Heisenberg spin chain for N-site PBC\n";
      cerr << "usage: heisenberg-energy <N>\n";
      return 1;
   }

   int const N = boost::lexical_cast<int>(argv[1]);

   int Niter = 100;  // guess - this needs to be bigger for large chains!

   int r = N/2;   // for the groundstate, we have N/2 spinons

   std::vector<double> z(r, 0.0);
   std::vector<double> zn(r, 0.0);
   
   for (int i = 0; i < Niter; ++i)
   {
      Iterate(N, zn, z);
      Iterate(N, z, zn);
   }

   // energy
   double Energy = 0;
   for (int i = 0; i < r; ++i)
   {
      //      cout << i << ' ' << Ir(N,i) << ' ' << z[i] << '\n';
      Energy += epsilon(z[i]);
   }

   cout.precision(16);
   cout << "Energy is " << ((Energy/N)+0.25) << endl;
}