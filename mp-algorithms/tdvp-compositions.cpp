// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tdvp-compositions.cpp
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "tdvp.h"
#include "tebd.h"

LTSDecomposition
SymmetricDecomposition(int Order, std::string Description,
                       std::vector<double> a,
                       std::vector<double> b)
{
   if (a.size() == b.size())
   {
      a.push_back(0.5 - std::accumulate(a.begin(), a.end(), 0.0));
      b.push_back(1.0 - 2.0 * std::accumulate(b.begin(), b.end(), 0.0));

      int asz = a.size();
      for (int i = asz-1; i >=0; --i)
         a.push_back(a[i]);

      int bsz = b.size();
      for (int i = bsz-2; i >= 0; --i)
         b.push_back(b[i]);
   }
   else if (a.size() == b.size()+1)
   {
      a.push_back(1.0 - 2.0 * std::accumulate(a.begin(), a.end(), 0.0));
      b.push_back(0.5 - std::accumulate(b.begin(), b.end(), 0.0));

      int asz = a.size();
      for (int i = asz-2; i >=0; --i)
         a.push_back(a[i]);

      int bsz = b.size();
      for (int i = bsz-1; i >= 0; --i)
         b.push_back(b[i]);
   }
   else
   {
      PANIC("Invalid SymmetricDecomposition");
   }

   return LTSDecomposition(Order, Description, a, b);
}

// The LeapfrogDecomposition assumes that the number of terms (m-1)/2 is odd, otherwise we need to repeat
// the final term again. We only need to supply n-1 weights, since the final weight is always
// 1 - 2 * sum(w)
// This is a product of leapfrog terms
// U(w) = e^{0.5*w*A} e^{w*B} r^{0.5*w*A}
// and for $n$ terms the expansion is U(w_1) U(w_2) ... U(w_{n-1}) U(w_n) U(w_{n-1}) ... U(w_1) ((m-1)/2 odd)
// or U(w_1) U(w_2) ... U(w_{n-1}) U(w_n) U(w_n) U(w_{n-1}) ... U(w_1) ((m-1)/2 even)
// We only need to supply (n-1) weights, as the final weight is 1 - 2*sum_{i<n}(w_i)
LTSDecomposition
LeapfrogDecompositionOdd(int Order, std::string Description, std::vector<double> w)
{
   std::vector<double> a;
   std::vector<double> b;
   //w.push_back(1.0 - 2.0 * std::accumulate(w.begin(), w.end(), 0.0));
   double aa = 0;
   for (auto ww : w)
   {
      a.push_back(0.5*(aa+ww));
      b.push_back(ww);
      aa = ww;
   }
   return SymmetricDecomposition(Order, Description, a, b);
}

Composition
ToComposition(std::string Description, LTSDecomposition d)
{
   std::vector<double> a = d.a();
   std::vector<double> b = d.b();

   auto A = a.begin();
   auto B = b.begin();

   std::vector<double> Gamma;

   Gamma.push_back(*A);
   Gamma.push_back(*B - Gamma.back());
   ++A, ++B;

   while (A != a.end()-1)
   {
      Gamma.push_back(*A - Gamma.back());
      Gamma.push_back(*B - Gamma.back());
      ++A, ++B;
   }

   if (Description.empty())
      Description = d.description();

   return Composition(d.order(), Description, Gamma);
}
