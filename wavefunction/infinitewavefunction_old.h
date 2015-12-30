// -*- C++ -*-

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
   return InfiniteWavefunctionLeft(Psi.PsiLinear, Psi.C_right, Psi.QShift);
}

#endif
