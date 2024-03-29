// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mps/packunpack.h
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
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
//
// pack & unpack functions for MatrixOperator.  These convert to/from
// an inner-product conserving vector representation.
// Preserving the inner product means we need an SU(2) coefficient.
//

#if !defined(MPTOOLKIT_MPS_PACKUNPACK_H)
#define MPTOOLKIT_MPS_PACKUNPACK_H

#include "state_component.h"  // for MatrixOperator definition
#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include <complex>

class PackMatrixOperator
{
   public:
      typedef std::complex<double> value_type;

      PackMatrixOperator(VectorBasis const& Basis1,
                         VectorBasis const& Basis2,
                         QuantumNumbers::QuantumNumber const& q);

                         PackMatrixOperator(VectorBasis const& Basis1,
                                            VectorBasis const& Basis2);

      PackMatrixOperator(MatrixOperator const& m);

      value_type* pack(MatrixOperator const& m, value_type* Iter) const;

      MatrixOperator unpack(value_type const* Iter) const;

      std::size_t size() const { return Size_; }

   private:
      void Initialize(VectorBasis const& Basis1,
                      VectorBasis const& Basis2,
                      QuantumNumbers::QuantumNumber const& q_);

      // OffsetRecType contains the information for an individual dense block
      struct OffsetRecType
      {
         int r, c;          // row and column of the MatrixOperator
         int Offset;        // offset into the packed array of this section
         unsigned Size;     // size of this block
         OffsetRecType() {}
         OffsetRecType(int r_, int c_, int Offset_, int Size_)
            : r(r_), c(c_), Offset(Offset_), Size(Size_) {}
      };

      typedef  std::vector<OffsetRecType> OffsetArrayType;

      VectorBasis B1_, B2_;
      QuantumNumbers::QuantumNumber q_;
      int Size_;  // linear size
      std::vector<OffsetRecType> OffsetArray_;   // array of blocks
      LinearAlgebra::Matrix<int> OffsetMatrix_;  // offset of each valid block
};

class PackStateComponent
{
   public:
      typedef std::complex<double> value_type;

      PackStateComponent(BasisList const& LocalBasis, VectorBasis const& Basis1,
                         VectorBasis const& Basis2);

      PackStateComponent(StateComponent const& m);

      value_type* pack(StateComponent const& m, value_type* Iter) const;

      StateComponent unpack(value_type const* Iter) const;

      std::size_t size() const { return Size_; }

   private:
      void Initialize(BasisList const& LocalBasis, VectorBasis const& Basis1,
                      VectorBasis const& Basis2);

      // OffsetRecType contains the information for an individual dense block
      struct OffsetRecType
      {
         int q;             // index into the StateComponent
         int r, c;          // row and column of the StateComponent
         int Offset;        // offset into the packed array of this section
         unsigned Size;     // size of this block
         OffsetRecType() {}
         OffsetRecType(int q_, int r_, int c_, int Offset_, int Size_)
            : q(q_), r(r_), c(c_), Offset(Offset_), Size(Size_) {}
      };

      typedef  std::vector<OffsetRecType> OffsetArrayType;

      BasisList B_;
      VectorBasis B1_, B2_;
      int Size_;  // linear size
      std::vector<OffsetRecType> OffsetArray_;   // array of blocks
      typedef std::vector<LinearAlgebra::Matrix<int>> OffsetMatrixType;
      OffsetMatrixType OffsetMatrix_;  // offset of each valid block
};

// some helper functors

template <typename F>
struct ApplyToPackedOperator
{
   ApplyToPackedOperator(F const& f, VectorBasis const& B1, VectorBasis const& B2,
			 QuantumNumbers::QuantumNumber const& q)
      : f_(f),
	Pack(B1, B2, q)
   {}

   ApplyToPackedOperator(F const& f, MatrixOperator const& M)
      : f_(f),
	Pack(M.Basis1(), M.Basis2(), M.TransformsAs())
   {}

   int pack_size() const { return Pack.size(); }

   void pack(MatrixOperator const& x, std::complex<double>* Out) const
   {
      Pack.pack(x, Out);
   }

   MatrixOperator unpack(std::complex<double> const* In) const
   {
      return Pack.unpack(In);
   }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = f_(x);
      Pack.pack(x, Out);
   }

   PackMatrixOperator Pack;
   F f_;
};

template <typename F>
struct ApplyToPackedStateComponent
{
   ApplyToPackedStateComponent(F const& f, BasisList const& B,
			       VectorBasis const& B1, VectorBasis const& B2)
      : f_(f),
	Pack(B, B1, B2)
   {}

   ApplyToPackedStateComponent(F const& f, StateComponent const& M)
      : f_(f),
	Pack(M.LocalBasis(), M.Basis1(), M.Basis2())
   {}

   int pack_size() const { return Pack.size(); }

   void pack(StateComponent const& x, std::complex<double>* Out) const
   {
      Pack.pack(x, Out);
   }

   StateComponent unpack(std::complex<double> const* In) const
   {
      return Pack.unpack(In);
   }
   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      StateComponent x = Pack.unpack(In);
      x = f_(x);
      Pack.pack(x, Out);
   }

   PackStateComponent Pack;
   F f_;
};

// construct the full matrix representation of some superoperator
// given by the functor F : MatrixOperator -> MatrixOperator
template <typename F>
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, MatrixOperator const& Init);

template <typename F>
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   TRACE(Size)(Size*Size);

   LinearAlgebra::Matrix<std::complex<double> > Out(Size, Size);

   std::vector<std::complex<double> > L(Size), R(Size);
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
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, StateComponent const& Init);

template <typename F>
LinearAlgebra::Matrix<std::complex<double> >
ConstructSuperOperator(F f, StateComponent const& Init)
{
   PackStateComponent Pack(Init);

   std::size_t Size = Pack.size();

   LinearAlgebra::Matrix<std::complex<double> > Out(Size, Size);

   std::vector<std::complex<double> > L(Size), R(Size);
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
LinearAlgebra::SparseMatrix<std::complex<double> >
ConstructSuperOperatorSparse(F f, MatrixOperator const& Init)
{
   PackMatrixOperator Pack(Init);

   std::size_t Size = Pack.size();

   LinearAlgebra::SparseMatrix<std::complex<double> > Out(Size, Size);

   std::vector<std::complex<double> > L(Size), R(Size);
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
         if (norm_frob(Out(j,i)) != 0)
            Out(j,i) = L[j];
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

   std::vector<std::complex<double> > L(Size), R(Size);
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
