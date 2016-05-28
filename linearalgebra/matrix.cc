// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrix.cc
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

namespace LinearAlgebra
{

//
// MatrixDimensions
//

#if defined(USE_PSTREAM)
template <int Format>
inline
PStream::opstreambuf<Format>& 
operator<<(PStream::opstreambuf<Format>& out, MatrixDimensions const& d)
{
   typename PStream::opstreambuf<Format>::size_type sz1 = d.size1, sz2 = d.size2;
   typename PStream::opstreambuf<Format>::difference_type s1 = d.stride1, s2 = d.stride2;
   return out << sz1 << sz2 << s1 << s2;
}

template <int Format>
inline
PStream::ipstreambuf<Format>& 
operator>>(PStream::ipstreambuf<Format>& in, MatrixDimensions& d)
{
   typename PStream::opstreambuf<Format>::size_type sz1, sz2;
   typename PStream::opstreambuf<Format>::difference_type s1, s2;
   in >> sz1 >> sz2 >> s1 >> s2;
   d.size1 = sz1;
   d.size2 = sz2;
   d.stride1 = s1;
   d.stride2 = s2;
   return in;
}
#endif

//
// MatrixRef
//

template <class Scalar, typename Orientation, typename Derived>
void 
MatrixRef<Scalar, Orientation, Derived>::resize(size_type NewRows, size_type NewCols)
{
   // Early return in the case that the data block is not shared and we don't need a resize.
   if (NewRows == this->size1() && NewCols == this->size2() 
       && this->stride2() == 1 && !Block.is_shared()) return;

   // this->cow(); // FIXME: pretty sure this is not needed???

   block_type NewData(NewRows*NewCols, NoHeader(), 
		      MatrixDimensions(NewRows, NewCols, Orientation()));
   Block = NewData;
}

template <class Scalar, typename Orientation, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixRef<Scalar, Orientation, Derived>&>::type
MatrixRef<Scalar, Orientation, Derived>::operator=(U const& x)
{
   // TODO: need a better temporary type
   Matrix<Scalar, Orientation> Temp(x);
   assign(*this, Temp);
   return *this;
}

template <class Scalar, typename Orientation, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixRef<Scalar, Orientation, Derived>&>::type
MatrixRef<Scalar, Orientation, Derived>::operator=(NoAliasProxy<U> const& x)
{
   assign(*this, x.value());
   return *this;
}

//
// Matrix
//

template <class Scalar, typename Orientation>
inline
Matrix<Scalar, Orientation>::Matrix(size_type Rows, size_type Cols)
  : MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >(Rows, Cols)
{
}

template <class Scalar, typename Orientation>
inline
Matrix<Scalar, Orientation>::Matrix(size_type Rows, size_type Cols, Scalar const& fill_)
  : MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >(Rows, Cols, fill_)
{
}

template <class Scalar, typename Orientation>
template <typename U>
inline
Matrix<Scalar, Orientation>::Matrix(U const& x, 
		       typename boost::enable_if<is_matrix<U> >::type*)
  : MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >(Size1<U>()(x), Size2<U>()(x))
{
   assign(*this, x);
}

template <class Scalar, typename Orientation>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, Matrix<Scalar, Orientation>&>::type
Matrix<Scalar, Orientation>::operator=(U const& x)
{
   // We must use a temporary here, as the right hand side could
   // be aliased without increasing the reference count
   // TODO: need a better temporary matrix type
   Matrix<Scalar, Orientation> Temp(x);
   this->resize(Temp.size1(), Temp.size2());
   assign(*this, Temp);
   return *this;
}

template <class Scalar, typename Orientation>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, Matrix<Scalar, Orientation>&>::type
Matrix<Scalar, Orientation>::operator=(NoAliasProxy<U> const& x)
{
   this->resize(Size<typename NoAliasProxy<U>::value_type>()(x.value()));
   assign(*this, x.value());
   return *this;
}

template <class Scalar, typename Orientation>
inline
Matrix<Scalar, Orientation>&
Matrix<Scalar, Orientation>::operator=(Matrix<Scalar, Orientation> const& V)
{
   this->block() = V.block().copy();
   return *this;
}

} // namespace LinearAlgebra
