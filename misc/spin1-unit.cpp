// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/spin1-unit.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "siteoperator/siteoperator.h"
#include "quantumnumbers/su2.h"
#include "common/math_const.h"
#include <map>

using namespace QuantumNumbers;
using namespace LinearAlgebra;

typedef std::complex<double> complex;

double factorial(double x)
{
   if (x == 1) return 1;
   else return x * factorial(x-1);
}

struct RedOperator
{
   RedOperator() {}

   RedOperator(SiteOperator S)
   {
      Op[S.TransformsAs()] = S;
   }

   std::map<QuantumNumber, SiteOperator> Op;
};

RedOperator operator*(RedOperator const& r, SiteOperator const& s)
{
   RedOperator Result;
   for (std::map<QuantumNumber, SiteOperator>::const_iterator
           I = r.Op.begin(); I != r.Op.end(); ++I)
   {
      QuantumNumberList Q = transform_targets(I->first, s.TransformsAs());
      for (std::size_t i = 0; i < Q.size(); ++i)
      {
         SiteOperator x = prod(I->second, s, Q[i]);
         if (norm_frob_sq(x) > 1E-12)
         {
            if (Result.Op.count(Q[i]))
               Result.Op[Q[i]] += x;
            else
               Result.Op[Q[i]] = x;
         }
      }
   }
   return Result;
}

RedOperator operator*(complex x, RedOperator const& r)
{
   RedOperator Result = r;
   for (std::map<QuantumNumber, SiteOperator>::iterator
           I = Result.Op.begin(); I != Result.Op.end(); ++I)
   {
      I->second *= x;
   }
   return Result;
}

RedOperator operator+(RedOperator const& x, RedOperator const& r)
{
   RedOperator Result = x;
   for (std::map<QuantumNumber, SiteOperator>::const_iterator
           I = r.Op.begin(); I != r.Op.end(); ++I)
   {
      if (Result.Op.count(I->first))
         Result.Op[I->first] += I->second;
      else
         Result.Op[I->first] = I->second;
   }
   return Result;
}

RedOperator operator*(RedOperator const& r, RedOperator const& s)
{
   RedOperator Result;
   for (std::map<QuantumNumber, SiteOperator>::const_iterator
           I = s.Op.begin(); I != s.Op.end(); ++I)
   {
      Result = Result + r*I->second;
   }
   return Result;
}

std::ostream& operator<<(std::ostream& out, RedOperator const& r)
{
   for (std::map<QuantumNumber, SiteOperator>::const_iterator
           I = r.Op.begin(); I != r.Op.end(); ++I)
   {
      out << "Component transforms as " << I->first << ":\n";
      out << I->second << '\n';
   }
   return out;
}

RedOperator exp(SiteOperator I, SiteOperator const& S)
{
   RedOperator e(I);
   RedOperator si(I);
   complex Factor = complex(0,math_const::pi);
   for (int i = 1; i < 50; ++i)
   {
      si = (1.0 / double(i))*si*S;
      e = e + si;
   }
   return e;
}

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);

   half_int Spin = 1;

   SiteBasis Basis(Symmetry);
   Basis.push_back("1", QN(1));

   SiteOperator S = SiteOperator(Basis, QN(1), SiteOperator::Bosonic);
   SiteOperator I = SiteOperator(Basis, QN(0), SiteOperator::Bosonic);
   S("1", "1") = sqrt(Spin * (Spin+1));
   I("1", "1") = 1.0;
   //S("1", "1") = 1; // sqrt(Spin * (Spin+1));

   SiteOperator Q = prod(S, S, QN(2));

#if 0
   TRACE(S);
   TRACE(Q);
   TRACE(prod(S, S, QN(0)));
   TRACE(prod(S, S, QN(1)));
   TRACE(prod(S, S, QN(2)));

   TRACE(prod(S, Q, QN(1)));
   TRACE(prod(S, Q, QN(2)));
   TRACE(prod(S, Q, QN(3)));

   TRACE(prod(Q, Q, QN(0)));
   TRACE(prod(Q, Q, QN(1)));
   TRACE(prod(Q, Q, QN(2)));
   TRACE(prod(Q, Q, QN(3)));
   TRACE(prod(Q, Q, QN(4)));
#endif

   std::cout.precision(12);

   for (double theta = -4; theta < 4; theta += 0.01)
   {
      std::cout << theta << '\n' << exp(I, theta*S) << "\n\n";
   }
}
