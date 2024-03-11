// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/mpsvd.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPSVD_H_DHCJKDSHTUY7893Y789YT578YGTE7R8OFDSHL)
#define MPSVD_H_DHCJKDSHTUY7893Y789YT578YGTE7R8OFDSHL

#include "density.h"

inline
void SingularValueDecomposition(SimpleOperator const& M,
                                SimpleOperator& U,
                                LinearAlgebra::Vector<double>& D,
                                SimpleOperator& Vt)
{
   CHECK(is_scalar(M.TransformsAs()));
   LinearBasis<BasisList> b1(M.Basis1()), b2(M.Basis2());
   CHECK_EQUAL(b1.size(), b2.size());

   typedef LinearAlgebra::Matrix<std::complex<double> > MatType;

   std::vector<MatType> MatList(b1.size());

   // The linear basis is ordered, so here i is iterating over quantum numbers
   // in both b1 and b2
   // Also construct the middle basis
   BasisList MiddleBasis(M.GetSymmetryList());
   for (std::size_t i = 0; i < b1.size(); ++i)
   {
      MatList[i] = MatType(b1.dim(i), b2.dim(i), 0.0);
      for (int n = 0; n < std::min(b1.dim(i), b2.dim(i)); ++n)
      {
         MiddleBasis.push_back(b1[i]);
      }
   }

   U = SimpleOperator(M.Basis1(), MiddleBasis, M.TransformsAs());
   D = LinearAlgebra::Vector<double>(MiddleBasis.size());
   Vt = SimpleOperator(MiddleBasis, M.Basis2(), M.TransformsAs());

   for (const_iterator<SimpleOperator>::type I = iterate(M); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         int si, ni, sj, nj;
         std::tie(si, ni) = b1.Lookup(J.index1());
         std::tie(sj, nj) = b2.Lookup(J.index2());
         CHECK_EQUAL(si,sj);
         MatList[si](ni,nj) = *J;
      }
   }

   int n = 0; // index of singular values, also index into the MiddleBasis
   for (std::size_t i = 0; i < b1.size(); ++i)
   {
      MatType Um, Vtm;
      LinearAlgebra::Vector<double> Dm;
      LinearAlgebra::SingularValueDecomposition(MatList[i], Um, Dm, Vtm);
      for (std::size_t k = 0; k < size(Dm); ++k)
      {
         for (std::size_t j = 0; j < size1(Um); ++j)
         {
            int r = b1.ReverseLookup(i, j);
            U(r,n) = Um(j,k);
         }
         for (std::size_t l = 0; l < size2(Vtm); ++l)
         {
            int r = b2.ReverseLookup(i, l);
            Vt(n,r) = Vtm(k,l);
         }
         D[n] = Dm[k];
         ++n;  // next singular value
      }
   }
}

#endif
