// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/infinitewavefunction_old.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTION_OLD_H)

#include "infinitewavefunctionleft.h"
#include "pstream/pstream.h"


struct InfiniteWavefunctionOld
{
      MatrixOperator C_old;
      QuantumNumbers::QuantumNumber QShift;
      LinearWavefunction PsiLinear;
      MatrixOperator C_right;
      AttributeList Attr;
};

inline
PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteWavefunctionOld& Psi)
{
   in >> Psi.C_old;
   in >> Psi.QShift;
   in >> Psi.PsiLinear;
   in >> Psi.C_right;
   in >> Psi.Attr;

   return in;
}

inline
PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionOld const& Psi)
{
   out << Psi.C_old;
   out << Psi.QShift;
   out << Psi.PsiLinear;
   out << Psi.C_right;
   out << Psi.Attr;
   return out;
}

inline
InfiniteWavefunctionLeft Make(InfiniteWavefunctionOld const& Psi)
{
   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(Psi.PsiLinear, Psi.C_right, Psi.QShift);
}

#endif
