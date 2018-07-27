// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix_view.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_VIEW_H)
#define MPTOOLKIT_BLAS_MATRIX_VIEW_H

#include "matrixref.h"

namespace blas
{

template <typename ValueType, typename Tag>
class matrix_range_view : public NormalMatrixProxy<ValueType, matrix_range_view<ValueType, Tag>, Tag>
{
   public:
      using derived_type       = matrix_range_view;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<ValueType>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;
      using remove_proxy_t     = Matrix<ValueType, Tag>;

      matrix_range_view(int Rows_, int Cols_, int LeadingDimension_, storage_type Ptr_)
         : Rows(Rows_), Cols(Cols_), LeadingDimension(LeadingDimension_), Ptr(Ptr_) {}

      matrix_range_view(matrix_range_view const&) = default;
      matrix_range_view(matrix_range_view&&) = default;

      ~matrix_range_view() = default;

      template <typename U>
      matrix_range_view& operator=(MatrixRef<ValueType, U, Tag> const& E) &&
      {
         assign(static_cast<matrix_range_view&&>(*this), E.as_derived());
         return *this;
      }

      template <typename U>
      matrix_range_view& operator+=(MatrixRef<ValueType, U, Tag> const& E) &&
      {
         add(static_cast<matrix_range_view&&>(*this), E.as_derived());
         return *this;
      }

      template <typename U>
      matrix_range_view& operator-=(MatrixRef<ValueType, U, Tag> const& E) &&
      {
         subtract(static_cast<matrix_range_view&&>(*this), E.as_derived());
         return *this;
      }

      int rows() const { return Rows; }
      int cols() const { return Cols; }
      int leading_dimension() const { return LeadingDimension; }

      storage_type storage() & { return Ptr; }
      const_storage_type storage() const& { return Ptr; }

   private:
      int Rows;
      int Cols;
      int LeadingDimension;
      storage_type Ptr;
};

} // namespace blas

#endif
