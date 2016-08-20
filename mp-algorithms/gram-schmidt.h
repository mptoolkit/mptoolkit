// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/gram-schmidt.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include <complex>
#include <vector>
#include "linearalgebra/vector.h"

namespace LinearSolvers
{

using namespace LinearAlgebra;
typedef std::complex<double> complex;

template <typename VectorType>
bool GramSchmidtAppend(std::vector<VectorType>& Basis,
                       VectorType NewVec,
                       double Ortho = 2.0, double ParallelThreshold = 0.0)
{
   int BasisSize = Basis.size();
   Basis.push_back(NewVec);

   double OriginalNorm = norm_frob_sq(NewVec);
   double Norm2 = OriginalNorm;
   bool Converged = false;
   while (!Converged)
   {
      double MaxOverlap = 0.0;
      for (int i = 0; i < BasisSize; ++i)
      {
         complex Overlap = inner_prod(Basis[i], Basis.back());
         MaxOverlap = std::max(norm_frob_sq(Overlap), MaxOverlap);
         Basis.back() -= Overlap * Basis[i];
      }
      double NewNorm2 = norm_frob_sq(Basis.back());

      if (NewNorm2 / OriginalNorm <= ParallelThreshold)  // parallel - cannot add the vector.
      {
         Basis.pop_back();
         return false;
      }

      Converged = (MaxOverlap <= Ortho * sqrt(Norm2));
      Norm2 = NewNorm2;
   }
   Basis.back() *= (1.0 / sqrt(Norm2));  // normalize
   return true;
}

} // namespace LinearSolvers
