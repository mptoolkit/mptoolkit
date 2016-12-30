// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/packunpack.h
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

#if !defined(PACKUNPACK_H_ASCHUIY243Y478Y648E89)
#define PACKUNPACK_H_ASCHUIY243Y478Y648E89

#include "mpstate.h"  // for MatrixOperator definition
#include "linearalgebra/vector.h"
#include <complex>

class PackMatrixOperator
{
   public:
      typedef std::complex<double> value_type;

      PackMatrixOperator(VectorBasis const& Basis1,
                         VectorBasis const& Basis2,
                         QuantumNumbers::QuantumNumber const& q);

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

#endif
