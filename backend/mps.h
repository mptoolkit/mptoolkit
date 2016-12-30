// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// backend/mps.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// An optimized representation of an MPS for DMRG and evolution algorithms
// This offers a limited range of operations - basically, designed for
// high performance DMRG.
//
// Needed operations:
// Construction via E matrices and F matrices
// inner_prod
// norm_frob
// matrix-vector multiply
//
// We'll also need some other operations, eg SVD
//
// As a first step, perhaps do construction from an existing
// StateComponent or MatrixOperator.
//
// The same structure is probably useable in both centre-matrix
// and centre-site formulations.


template <typename Real>
class StateDescription
{
   public:

   private:



};



struct InnerProdElement;
{
   int Offset;
   int Size;
   double Prefactor;
};

typedef std::vector<InnerProdElement> InnerProdDescription;

double inner_prod(InnerProdDescription const& Desc,
                  double const* v1, double const* v2)
{
   double Result = 0.0;
   unsigned sz = Desc.size();
   for (unsigned i = 0; i < sz; ++i)
   {
      Result += Desc[i].Prefactor
         * inner_prod(v1+Desc[i].Offset,
                      v2+Desc[i].Offset,
                      Desc[i].Size);
   }
   return Result;
}


double norm_frob(InnerProdDescription const& Desc,
                 double const* v)
{
   double Result = 0.0;
   unsigned sz = Desc.size();
   for (unsigned i = 0; i < sz; ++i)
   {
      Result += Desc[i].Prefactor
         * norm_frob(v+Desc[i].Offset, Desc[i].Size);
   }
   return Result;
}

// typedef for the integer type for describing offsets relative
// to a memory block
typedef ptrdiff_t offset_t;

// A description of a block of an MPS A^{Local}_{Row,Column}
// Local is an index into the LocalBasis
// Row is an index into Basis1
// Column is an index into Basis2
struct ComponentDescription
{
   int Local;
   int Row;
   int Column;

   // Pointer into the storage location
   offset_t Offset;
};


struct StateComponent
{

   BasisList LocalBasis_;
   VectorBasis Basis1_;
   VectorBasis Basis2_;
   std::vector<ComponentDescription> Components_;
   double* Storage_;

   InnerProdDescription InnerProd_;
};
