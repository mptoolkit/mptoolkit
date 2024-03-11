// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/vectorfwd.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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

/*
  vectorfwd.h

  Forward declarations so we can define the default value_types of expressions
  in vectorinterface.h

  Created 2005-01-12 Ian McCulloch
*/

#if !defined(VECTORFWD_H_CJKHUTYIUTYT897YT9874O)
#define VECTORFWD_H_CJKHUTYIUTYT897YT9874O

namespace LinearAlgebra
{

// dense vectors default to Vector
template <typename Scalar>
class Vector;

// sparse vectors default to Vector
template <typename T>
class MapVector;

} // namespace LinearAlgebra

#endif
