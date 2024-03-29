// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixblock.cc
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

namespace LinearAlgebra
{

template <class Scalar, typename Orientation, typername Header, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixBlockRef<Scalar, Orientation, Header, Derived>&>::type
MatrixBlockRef<Scalar, Orientation, Header, Derived>::operator=(U const& x)
{
   // TODO: need a better temporary type
   Matrix<Scalar, Orientation> Temp(x);
   assign(*this, Temp);
   return *this;
}

template <class Scalar, typename Orientation, typename Header, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixBlockRef<Scalar, Orientation, Header, Derived>&>::type
MatrixBlockRef<Scalar, Orientation, Header, Derived>::operator=(NoAliasProxy<U> const& x)
{
   assign(*this, x.value());
   return *this;
}

} // namespace LinearAlgebra
