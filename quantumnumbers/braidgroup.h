// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/braidgroup.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_QUANTUMNUMBERS_BRAIDGROUP_H)
#define MPTOOLKIT_QUANTUMNUMBERS_BRAIDGROUP_H

class BraidGroup
{
};

PStream::opstream& operator<<(PStream::opstream& out, BraidGroup const& x);
PStream::opstream& operator>>(PStream::opstream& out, BraidGroup& x);

#endif
