// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/braidrep.h
//
// Copyright (C) 2016-2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// BreadRep denotes a representation of a braid group.  There is one
// representation for each symmetry.  For example, with U(1) symmetries
// we have three classes of representations,
// bosonic, fermionic, or anyonic (with some phase theta).
// If we had a SymmetryList of "N:U(1)", then we might define
// a BraidGroup as eg "N:Bosonic", or "N:Fermionic", or "N:AbelianAnyon($theta)".
// A BraidRep then defines a specific representation of that braid group.

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_BRAIDREP_H)
#define MPTOOLKIT_QUANTUMNUMBERS_BRAIDREP_H

class BraidRep
{
   public:

   private:
};

PStream::opstream& operator<<(PStream::opstream& out, BraidRep const& x);
PStream::opstream& operator>>(PStream::opstream& out, BraidRep& x);

#endif
