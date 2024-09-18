// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrix_temp.h
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

/*
  matrix_temp.h

  A 'temporary' matrix using stackallocator.h

  Created 2004-05-14 Ian McCulloch
*/

#if !defined(MATRIX_TEMP_H_HSFUJRHT7843Y7FHOL894Y39HWH)
#define MATRIX_TEMP_H_HSFUJRHT7843Y7FHOL894Y39HWH

#include "matrixtypes.h"
#include "common/stackallocator.h"

namespace LinearAlgebra
{

//
// TempMatrixRef
//

template <typename Scalar>
class TempMatrixRef;

template <class Scalar>
struct MatrixTraits<TempMatrixRef<Scalar> >
{
   typedef MatrixConstSlice<Scalar, TempMatrixRef<Scalar> >          rvalue_type;
   typedef MatrixOuterIterator<Scalar*, Slice, Range>                iterator1;
   typedef MatrixOuterIterator<Scalar*, Range, Slice>                iterator2;
   typedef MatrixOuterIterator<Scalar const*, Slice, Range>          const_iterator1;
   typedef MatrixOuterIterator<Scalar const*, Range, Slice>          const_iterator2;
   typedef TempMatrixRef<Scalar>                                     container_reference;
  //   typedef TempMatrixRef<Scalar> const&                              container_reference;
   typedef TempMatrixRef<Scalar> const&                              expression_type;
   typedef TempMatrixRef<Scalar>                                     container_const_reference;
   typedef Scalar&                                                   result_type;
   typedef Scalar const&                                             const_result_type;

   typedef MatrixSection<Scalar, Range, Range, TempMatrixRef<Scalar> >      range1_type;
   typedef MatrixSection<Scalar, Slice, Range, TempMatrixRef<Scalar> >      slice1_type;
   typedef MatrixSection<Scalar, Index, Range, TempMatrixRef<Scalar> >      index1_type;
   typedef MatrixSection<Scalar, Range, Range, TempMatrixRef<Scalar> >      range2_type;
   typedef MatrixSection<Scalar, Range, Slice, TempMatrixRef<Scalar> >      slice2_type;
   typedef MatrixSection<Scalar, Range, Index, TempMatrixRef<Scalar> >      index2_type;

   typedef MatrixConstSection<Scalar, Range, Range, TempMatrixRef<Scalar> > const_range1_type;
   typedef MatrixConstSection<Scalar, Slice, Range, TempMatrixRef<Scalar> > const_slice1_type;
   typedef MatrixConstSection<Scalar, Index, Range, TempMatrixRef<Scalar> > const_index1_type;
   typedef MatrixConstSection<Scalar, Range, Range, TempMatrixRef<Scalar> > const_range2_type;
   typedef MatrixConstSection<Scalar, Range, Slice, TempMatrixRef<Scalar> > const_slice2_type;
   typedef MatrixConstSection<Scalar, Range, Index, TempMatrixRef<Scalar> > const_index2_type;

   typedef MatrixTransposeProxy<Scalar, TempMatrixRef<Scalar> >            transpose_type;
   typedef MatrixConstTransposeProxy<Scalar, TempMatrixRef<Scalar> >       const_transpose_type;
};

template <typename Scalar>
class TempMatrixRef : public MatrixSlice<Scalar, TempMatrixRef<Scalar> >
{
   public:
      typedef MatrixSlice<Scalar, TempMatrixRef<Scalar> >   base_class;

      typedef MatrixTraits<TempMatrixRef<Scalar> >          traits_type;

      typedef typename traits_type::iterator1            iterator1;
      typedef typename traits_type::iterator2            iterator2;
      typedef typename traits_type::const_iterator1      const_iterator1;
      typedef typename traits_type::const_iterator2      const_iterator2;
      typedef typename traits_type::transpose_type       transpose_type;
      typedef typename traits_type::const_transpose_type const_transpose_type;

      TempMatrixRef(TempMatrixRef const& V) : Size1(V.Size1), Size2(V.Size2), Data(V.Data) {}

      template <class S2, class D2>
      TempMatrixRef& operator=(MatrixConstExpression<S2, D2> const& V)
      {
         this->assign_copy(V.as_derived());
         return *this;
      }

      TempMatrixRef& operator=(TempMatrixRef const& V)
      {
         this->assign(V.as_derived());
         return *this;
      }

      size_type size() const { return Size1 * Size2; }
      size_type size1() const { return Size1; }
      size_type size2() const { return Size2; }

      iterator1 begin1() { return iterator1(Data, Slice(0, Size1, Size2), 0, Range(0, Size2)); }
      iterator1 end1() { return iterator1(Data, Slice(0, Size1, Size2), Size1, Range(0, Size2)); }

      iterator2 begin2() { return iterator2(Data, Range(0, Size2), 0, Slice(0, Size1, Size2)); }
      iterator2 end2() { return iterator2(Data, Range(0, Size2), Size2, Slice(0, Size1, Size2)); }

      const_iterator1 begin1() const
        { return const_iterator1(Data, Slice(0, Size1, Size2), 0, Range(0, Size2)); }

      const_iterator1 end1() const
        { return const_iterator1(Data, Slice(0, Size1, Size2), Size1, Range(0, Size2)); }

      const_iterator2 begin2() const
        { return const_iterator2(Data, Range(0, Size2), 0, Slice(0, Size1, Size2)); }

      const_iterator2 end2() const
        { return const_iterator2(Data, Range(0, Size2), Size2, Slice(0, Size1, Size2)); }

      Scalar const& operator()(size_type i, size_type j) const
        { return Data[i * Size2 + j]; }

      Scalar& operator()(size_type i, size_type j)
        { return Data[i * Size2 + j]; }


      // members defined in MatrixStride
      difference_type stride1() const { return Size2; }
      difference_type stride2() const { return 1; }
      Scalar* data() { return Data; }
      Scalar const* data() const { return Data; }

   protected:
      TempMatrixRef(size_type Size1, size_type Size2);

      TempMatrixRef(size_type Size1, size_type Size2, Scalar const& Init);

      void deallocate();

   private:
      void init(size_type Size1_, size_type Size2_);

      size_type Size1, Size2;
      Scalar* Data;
};

#if 0
template <typename Scalar, typename Nest>
struct MatrixTranspose<TempMatrixRef<Scalar>, Nest>
   : public MatrixTransform<MatrixTransformProxy<TempMatrixRef<Scalar> >, Nest>
{
   MatrixTranspose(Nest n = Nest())
      : MatrixTransform<MatrixTransformProxy<TempMatrixRef<Scalar> >, Nest>(n) {}
};

template <typename Scalar, typename Nest>
struct MatrixConstTranspose<TempMatrixRef<Scalar>, Nest>
   : public MatrixTransform<MatrixConstTransformProxy<TempMatrixRef<Scalar> >, Nest>
{
   MatrixConstTranspose(Nest n = Nest())
      : MatrixConstTransform<MatrixTransformProxy<TempMatrixRef<Scalar> >, Nest>(n) {}
};
#endif

template <typename Scalar>
void TempMatrixRef<Scalar>::init(size_type Size1_, size_type Size2_)
{
   Size1 = Size1_;
   Size2 = Size2_;
   Data = static_cast<Scalar*>(StackAlloc::allocate(Size1 * Size2 * sizeof(Scalar)));
}

template <typename Scalar>
void TempMatrixRef<Scalar>::deallocate()
{
   StackAlloc::deallocate(Data, Size1*Size2*sizeof(Scalar));
}

template <typename Scalar>
TempMatrixRef<Scalar>::TempMatrixRef(size_type Size1, size_type Size2)
{
   this->init(Size1, Size2);
}

template <typename Scalar>
TempMatrixRef<Scalar>::TempMatrixRef(size_type Size1, size_type Size2, Scalar const& Init)
{
   this->init(Size1, Size2);
   this->fill(Init);
}

//
// TempMatrix
//

template <typename Scalar>
class TempMatrix : public TempMatrixRef<Scalar>
{
   public:
      typedef MatrixSlice<Scalar, TempMatrixRef<Scalar> >   base_class;

      typedef MatrixTraits<TempMatrixRef<Scalar> >          traits_type;

      typedef typename traits_type::iterator1            iterator1;
      typedef typename traits_type::iterator2            iterator2;
      typedef typename traits_type::const_iterator1      const_iterator1;
      typedef typename traits_type::const_iterator2      const_iterator2;
      typedef typename traits_type::transpose_type       transpose_type;
      typedef typename traits_type::const_transpose_type const_transpose_type;

      TempMatrix(size_type Size1, size_type Size2)
        : TempMatrixRef<Scalar>(Size1, Size2) {}

      TempMatrix(size_type Size1, size_type Size2, Scalar const& Init)
        : TempMatrixRef<Scalar>(Size1, Size2, Init) {}

      TempMatrix(TempMatrix const& V)
        : TempMatrixRef<Scalar>(V.size1(), V.size2()) { this->assign(V); }

      template <class S2, class D2>
      TempMatrix(GenericMatrix<S2, D2> const& V)
        : TempMatrixRef<Scalar>(V.size1(), V.size2()) { this->assign(V); }

      template <class S2, class D2>
      TempMatrix& operator=(GenericMatrix<S2, D2> const& V)
      {
         this->assign_copy(V.as_derived());
         return *this;
      }

      TempMatrix& operator=(TempMatrix const& V)
      {
         this->assign(V.as_derived());
         return *this;
      }

      ~TempMatrix() { this->deallocate(); }

};

} // namespace LinearAlgebra

#endif
