// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/denseoperator.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "denseoperator.h"

std::ostream& operator<<(std::ostream& out, DenseOperator const& x)
{
   out << "DenseOperator:\nBasis1() = " << x.Basis1()
       << "\nBasis2() = " << x.Basis2()
       << "data() =\n" << x.data();
   return out;
}
