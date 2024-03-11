// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/get-k.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER


#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <map>

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"

using namespace std;


std::complex<double> ComplexFromString(std::string const& s)
{
   char* End;
   double RealPart = strtod(s.c_str(), &End);
   double ImagPart = (End == &s[s.size()-1]) ? 0.0
      : strtod(End, &End);
   return std::complex<double>(RealPart,ImagPart);
}

typedef LinearAlgebra::Matrix<complex<double> > CMatrix;


LinearAlgebra::Vector<complex<double> >
Eigenvalues(CMatrix const& M)
{
   CMatrix L, R;
   LinearAlgebra::Vector<complex<double> > Result = Diagonalize(M, L, R);
   return Result;
#if 0


   if (size1(M) != 2)
      return LinearAlgebra::Vector<complex<double> >();

   complex<double> T = M(0,0) + M(1,1);

   complex<double> D = M(0,0)*M(1,1) - M(0,1)*M(1,0);

   complex<double> Delta = sqrt(0.25*T*T - D);

   LinearAlgebra::Vector<complex<double> > Result(2);
   Result[0] = (T/2.0) - D;
   Result[1] = (T/2.0) + D;
   return Result;
#endif
}

int main(int argc, char** argv)
{
   if (argc != 3)
   {
      cout << "expected: <densitymatrix> <transfermatrix>\n";
      return 1;
   }

   map<int, double> Density;
   ifstream In(argv[1]);
   int i, j;
   double x;
   while (In >> i >> j >> x)
   {
      if (i != j)
      {
         cerr << "error: density matrix should be diagonal.\n";
         return 1;
      }

      Density[i] = x;
   }

   // find near-degenerate eigenvalues

   double const Closeness = 0.1;

   typedef std::map<double, std::set<int> > DMGroupType;
   std::map<double, std::set<int> > DMGroup;

   for (map<int, double>::const_iterator I = Density.begin();
        I != Density.end(); ++I)
   {
      bool Found = false;
      // are we sufficently close in magnitude to an existing eigenvalue?
      for (DMGroupType::iterator J = DMGroup.begin(); J != DMGroup.end(); ++J)
      {
         if (I->second > (J->first * (1-Closeness)) && (I->second < (J->first * (1+Closeness))))
         {
            std::cout << "Eigenvalue " << I->second << " at row " << I->first
                      << " matches existing eigenvalue " << J->first << " represented at " <<
               (*J->second.begin()) << '\n';
            J->second.insert(I->first);
            Found = true;
            break;
         }
      }
      // otherwise, make a new group
      if (!Found)
      {
         std::cout << "New eigenvalue group " << I->second << " at row " << I->first << '\n';
         DMGroup[I->second] = std::set<int>();
         DMGroup[I->second].insert(I->first);
      }
   }

   // Make matrices
   vector<CMatrix> TM(DMGroup.size());
   vector<double> Eigenvalue;
   // initialize
   i = 0;
   for (DMGroupType::const_iterator I = DMGroup.begin(); I != DMGroup.end(); ++I)
   {
      TM[i] = CMatrix(I->second.size(), I->second.size(), 0.0);
      ++i;
      Eigenvalue.push_back(I->first);
   }

   // load transfer matrixes
   std::string Str;
   ifstream TransFile(argv[2]);
   while (TransFile >> i >> j >> Str)
   {
      complex<double> x = ComplexFromString(Str);
      // find the coordinates
      int n = 0;
      bool Set = false;
      for (DMGroupType::const_iterator I = DMGroup.begin(); I != DMGroup.end(); ++I)
      {
         int c1=-1, c2=-1;
         int Index = 0;
         for (std::set<int>::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
         {
            if (*J == i)
               c1 = Index;
            if (*J == j)
               c2 = Index;
            ++Index;
         }
         if (c1 != -1 && c2 != -1)
         {
            Set = true;
            cout << "Matrix element in block " << n << " at index " << c1 << "," << c2
                 << " is " << x << '\n';
            TM[n](c1,c2) = x;
            break;
         }
         ++n;
         if (Set)
            break;
      }
   }

   double const pi = 3.1415926535897;

   for (unsigned k = 0; k < TM.size(); ++k)
   {
      cout << "Eigenvalue group: " << Eigenvalue[k] << '\n';
      cout << TM[k] << "\n";
      cout << (TM[k] * herm(TM[k])) << "\n";
      LinearAlgebra::Vector<complex<double> > e = Eigenvalues(TM[k]);
      if (size(e) != 0)
      {
         cout << "momenta: ";
         for (unsigned n = 0; n < size(e); ++n)
         {
            if (n > 0)
               cout << ", ";
            cout << atan2(e[n].imag(), e[n].real())/pi;
         }
         cout << "\n";
      }
      cout << "\n";
   }

}
