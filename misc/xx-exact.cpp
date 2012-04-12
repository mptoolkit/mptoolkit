// -*- C++ -*- $Id$

/*
  Calculates the exact solution of the time evolution of the XX model
  starting from the initial state |up ... up up down down ... down>
  Link with GNU scientific library.
*/

#include <iostream>
#include <iomanip>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include <boost/lexical_cast.hpp>

using namespace std;

double XX_Sz_exact( double t, int n)
{ 
   double Sz=0;  int j0=n;
   if( j0<0 ) j0=-j0-1;
   for(int j=-j0; j<=j0; j++)
      Sz += -1./2 * pow( gsl_sf_bessel_Jn(j,t), 2);
   if(n<0) Sz=-Sz;
   return Sz;
}

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      cerr << "usage: xx-exact <L> <t>\n";
      return 1;
   }

   int L = boost::lexical_cast<int>(argv[1]);
   double t = boost::lexical_cast<double>(argv[2]);

   cout.precision(20);
   for (int n = 0; n < L; ++n)
   {
      cout << setw(25) << t << ' ' 
                << setw(25) <<  XX_Sz_exact(t, n-L/2) << endl;
   }
}
