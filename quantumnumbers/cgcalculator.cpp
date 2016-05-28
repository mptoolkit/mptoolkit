// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/cgcalculator.cpp
//
// Copyright (C) 2001-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

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
