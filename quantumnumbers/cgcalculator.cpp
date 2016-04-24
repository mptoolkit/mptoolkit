// -*- C++ -*- $Id$

/*
  simple Clebsch-Gordan coefficient calculator.

  Created 2001-08-19 Ian McCulloch
*/

#include <iostream>
#include <iomanip>
#include "coupling.h"
#include "common/bindpair.h"

using gmp::rational;
using namespace std;

int main()
{
  std::cout << "Clebsch-Gordan coefficient calculator.\n"
	    << "Enter spins as decimals, enter j1 < 0 to quit.\n\n";
  while (1)
  {
      half_int j1, j2, j, m1, m2, m;
      cout << "j1: "; cin >> j1;
      if (j1 < 0) return 0;

      cout << "j2: "; cin >> j2;
      cout << "j: "; cin >> j;
      
      for (half_int m = -j; m <= j; ++m)
      {
	 for (half_int m1 = -j1; m1 <= j1; ++m1)
         {
            for (half_int m2 = -j2; m2 <= j2; ++m2)
            {
	      if (m1+m2==m)
	      {
		 rational x, y;
		 std::tie(x,y) = ClebschGordanSquared(j1,m1,j2,m2,j,m);
		 rational a = x*x*y;
		 cout << "\nCG{ " << setw(4) << j1 << " , " << setw(4) << j2 << " , " << setw(4) << j << " }\n"
		      << "  { " << setw(4) << m1 << " , " << setw(4) << m2 << " , " << setw(4) << m << " } "
		      << " = ";
		 int s = gmp::sign(x);
		 if (s == -1) cout << "-sqrt(" << a.to_string() << ")\n\n";
		 else if (s == 0) cout << "0";
		 else if (s == 1) cout << "sqrt(" << a.to_string() << ")\n\n";
	      }
	    }
	 }
      }
  }
}
