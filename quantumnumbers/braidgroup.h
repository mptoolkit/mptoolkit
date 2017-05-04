// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/braidgroup.h
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
//
// A representation of the Braid group associated with a symmetry list.
// For example, with U(1) symmetries
// we have three classes of representations,
// bosonic, fermionic, or anyonic (with some phase theta).
// If we had a SymmetryList of "N:U(1)", then we might define
// a BraidGroup as eg "N:Bosonic", or "N:Fermionic", or "N:AbelianAnyon($theta)".

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_BRAIDGROUP_H)
#define MPTOOLKIT_QUANTUMNUMBERS_BRAIDGROUP_H

class BraidGroup
{
   public:

   private:
};

PStream::opstream& operator<<(PStream::opstream& out, BraidGroup const& x);
PStream::opstream& operator>>(PStream::opstream& out, BraidGroup& x);

#endif
