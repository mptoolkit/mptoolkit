// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/packunpack.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// pack & unpack functions for MatrixOperator.  These convert to/from
// an inner-product conserving vector representation.
// Preserving the inner product means we need an SU(2) coefficient.
//

#if !defined(MPTOOLKIT_MPS_PACKUNPACK_H)
#define MPTOOLKIT_MPS_PACKUNPACK_H

#include "state_component.h"  // for MatrixOperator definition
#include <complex>

// pack a MatrixOperator into an inner-product conserving linear vector
blas::Vector<complex>
pack(MatrixOperator const& m);

// unpack a linear vector into a MatrixOperator.  The operator
// must have the appropriate basis and quantum numbers set
void unpack(MatrixOperator& m, blas::Vector<complex> const& v);

blas::Vector<complex>
pack(StateComponent const& m);

void unpack(StateComponent& m, blas::Vector<complex> const& v);

// construct the full matrix representation of some superoperator
// given by the functor F : MatrixOperator -> MatrixOperator
template <typename F>
blas::Matrix<std::complex<double>>
ConstructSuperOperator(F f, MatrixOperator const& Init);

#if 0
// not updated yet to use the new pack/unpack API
template <typename F>
blas::Matrix<std::complex<double>>
ConstructSuperOperator(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   TRACE(Size)(Size*Size);

   blas::Matrix<std::complex<double> > Out(Size, Size);

   blas::Vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
         Out(j,i) = L[j];
      }
   }

   return Out;
}

// construct the full matrix representation of some superoperator
// given by the functor F : StateComponent -> StateComponent
template <typename F>
blas::Matrix<std::complex<double> >
ConstructSuperOperator(F f, StateComponent const& Init);

template <typename F>
blas::Matrix<std::complex<double> >
ConstructSuperOperator(F f, StateComponent const& Init)
{
   PackStateComponent Pack(Init);

   std::size_t Size = Pack.size();

   blas::Matrix<std::complex<double> > Out(Size, Size);

   blas::Vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      StateComponent M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
         Out(j,i) = L[j];
      }
   }

   return Out;
}


template <typename F>
blas::SparseMatrix<std::complex<double> >
ConstructSuperOperatorSparse(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   blas::SparseMatrix<std::complex<double> > Out(Size, Size);

   blas::Vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
         if (norm_frob(L[j]) != 0)
            Out.insert(j,i, L[j]);
      }
   }

   return Out;
}

// write the operator to a file, column major format
template <typename F>
std::streambuf*
ConstructSuperOperatorBinaryFile(F f, MatrixOperator const& Init, std::streambuf* Out)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   blas::Vector<std::complex<double> > L(Size), R(Size);
   for (unsigned i = 0; i < Size; ++i)
   {
      if (i > 0)
         R[i-1] = 0.0;

      R[i] = 1.0;

      MatrixOperator M = Pack.unpack(&R[0]);
      M = f(M);
      Pack.pack(M, &L[0]);
      for (unsigned j = 0; j < Size; ++j)
      {
         double r = L[j].real();
         double c = L[j].real();
         Out->sputn(reinterpret_cast<char const*>(&r), sizeof(r));
         Out->sputn(reinterpret_cast<char const*>(&c), sizeof(c));
      }
   }

   return Out;
}

#endif

#endif
