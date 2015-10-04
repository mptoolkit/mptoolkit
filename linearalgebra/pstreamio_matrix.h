/* -*- C++ -*- $Id$

  pstreamio_matrix.h

  I/O for matrices using the PStream framework.

  Created 2005-04-06 Ian McCulloch
*/ 

#if !defined(PSTREAMIO_MATRIX_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define PSTREAMIO_MATRIX_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "matrixoperations.h"

namespace LinearAlgebra
{

// The available formats for saving matrices.
// This isn't an enumeration because we need to control the representation
// on disk.
struct MatrixFormats
{
   typedef char type;
   static type const DenseRowMajor = 'm';
   static type const Coordinate = 'c';
   static type const Diagonal = 'd';
   static type const Fixed = 'f';
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Mat,
                             Concepts::DenseMatrix<MatV, RowMajor, MatI>>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Mat const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type M) const
   {
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Rows = size1(M), Cols = size2(M);
#if !defined(PSTREAMIO_OLD_FORMAT)
      out << MatrixFormats::DenseRowMajor;
#endif
      out << Rows << Cols;

      typename const_iterator<Mat>::type I = iterate(M);
      while (I)
      {
         typename const_inner_iterator<Mat>::type J = iterate(I);
         while (J)
         {
            out << *J;
            ++J;
         }
         ++I;
      }
      return out;
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Mat,
                             Concepts::DenseMatrix<MatV, ColMajor, MatI>>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Mat const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type M) const
   {
      return out << swap_sort_order(M);
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Mat,
                             Concepts::DiagonalMatrix<MatV, MatI>>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Mat const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type M) const
   {
      out << M.diagonal();
      return out;
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Mat,
                             Concepts::SparseMatrix<MatV, MatI>>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Mat const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type M) const
   {
#if defined(PSTREAMIO_OLD_FORMAT)
      PANIC("Cannot serialize sparse matrices with PSTREAMIO_OLD_FORMAT");
#endif
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Rows = size1(M), Cols = size2(M), Nnz = nnz(M);
      out << MatrixFormats::Coordinate;
      out << Rows << Cols << Nnz;

      typename const_iterator<Mat>::type I = iterate(M);
      while (I)
      {
         typename const_inner_iterator<Mat>::type J = iterate(I);
         while (J)
         {
            st r = J.index1(), c = J.index2();
            out << r << c;
            out << *J;
            ++J;
         }
         ++I;
      }
      return out;
   }
};



namespace Private
{

//
// implementation helpers
//

template <typename Stream, typename Mat>
void do_serialize_in_dense_row_major(Stream& in, Mat& M)
{
   typename iterator<Mat>::type I = iterate(M);
   while (I)
   {
      typename inner_iterator<Mat>::type J = iterate(I);
      while (J)
      {
         in >> *J;
         ++J;
      }
      ++I;
   }
}

template <typename Stream, typename Mat>
void do_serialize_in_dense_row_major(Stream& in, Mat const& M)
{
   do_serialize_in_dense_row_major(in, const_cast<Mat&>(M));
}

template <int Format, typename Mat>
void do_serialize_in_coordinate(PStream::ipstreambuf<Format>& in, Mat& M)
{
   typedef typename PStream::ipstreambuf<Format>::size_type st;
   typedef typename interface<Mat>::value_type value_t;
   st nnz;
   in >> nnz;
   for (st i = 0; i < nnz; ++i)
   {
      st r,c;
      in >> r >> c;
      set_element(M, r,c, in.template read<value_t>());
   }
}

template <int Format, typename Mat>
void do_serialize_in_coordinate(PStream::ipstreambuf<Format>& in, Mat const& M)
{
   do_serialize_in_coordinate(in, const_cast<Mat&>(M));
}

} // namespace Private

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Mat&,
                            Concepts::DenseMatrix<MatV, RowMajor, MatI>>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Mat& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type M) const
   {
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Rows, Cols;
      MatrixFormats::type t;
#if defined(PSTREAMIO_OLD_FORMAT)
      t = MatrixFormats::DenseRowMajor;
#else
      in >> t;
#endif
      if (t == MatrixFormats::DenseRowMajor)
      {
         in >> Rows >> Cols;
         try_resize(M, Rows, Cols);
         Private::do_serialize_in_dense_row_major(in, M);
      }
      else if (t == MatrixFormats::Coordinate)
      {
         in >> Rows >> Cols;
         try_resize(M, Rows, Cols);
         zero_all(M);
         Private::do_serialize_in_coordinate(in, M);
      }
      else
      {
         PANIC("Unsupported matrix format.")(t)(int(t));
      }
      return in;
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Mat&,
                            Concepts::DiagonalMatrix<MatV, MatI>>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Mat& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type M) const
   {
      in >> M.diagonal();
      return in;
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Mat&,
                            Concepts::DenseMatrix<MatV, ColMajor, MatI>>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Mat& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type M) const
   {
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Rows, Cols;
      MatrixFormats::type t;
#if defined(PSTREAMIO_OLD_FORMAT)
      t = MatrixFormats::DenseRowMajor;
#else
      in >> t;
#endif
      if (t == MatrixFormats::DenseRowMajor)
      {
         in >> Rows >> Cols;
         try_resize(M, Rows, Cols);
         Private::do_serialize_in_dense_row_major(in, swap_sort_order(M));
      }
      else if (t == MatrixFormats::Coordinate)
      {
         in >> Rows >> Cols;
         try_resize(M, Rows, Cols);
         zero_all(M);
         Private::do_serialize_in_coordinate(in, M);
      }
      else
      {
         PANIC("Unsupported matrix format.")(t)(int(t));
      }
      return in;
   }
};

template <int Format, typename Mat, typename MatV, typename MatI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Mat&,
                            Concepts::SparseMatrix<MatV, MatI>>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Mat& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type M) const
   {
#if defined(PSTREAMIO_OLD_FORMAT)
      PANIC("Cannot serialize sparse matrices with PSTREAMIO_OLD_FORMAT");
#endif
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Rows, Cols;
      MatrixFormats::type t;
      in >> t;
      if (t == MatrixFormats::Coordinate)
      {
         in >> Rows >> Cols;
         try_resize(M, Rows, Cols);
         zero_all(M);
         Private::do_serialize_in_coordinate(in, M);
      }
      else
      {
         PANIC("Unsupported matrix format for sparse matrix.")(t)(int(t));
      }
      return in;
   }
};

} // namespace LinearAlgebra

#endif
