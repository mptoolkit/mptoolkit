// -*- C++ -*- $Id$
//
// The wavefunction format from version 0.7.3, for backwards compatibility
//

#if !defined(MPWAVEFUNCTION_COMPAT_H_SDOJUIF3487R84789U89)
#define MPWAVEFUNCTION_COMPAT_H_SDOJUIF3487R84789U89

#include "linearwavefunction.h"

inline
LinearWavefunction LoadWavefunction_0_7_3(PStream::ipstream& in)
{
   typedef MPStateComponent ComponentType;
   typedef MatrixOperator OperatorType;
   // mpstate hasn't changed
   std::vector<ComponentType> LeftStack, RightStack;
   ComponentType LeftTop, RightTop;
   OperatorType CenterOp;
   AttributeList Attr;
   
   in >> LeftStack >> RightStack >> LeftTop >> RightTop >> CenterOp >> Attr;

   LinearWavefunction Psi(CenterOp.GetSymmetryList());
   Psi.Attributes() = Attr;
   Psi.push_back(RightTop);
   for (int i = RightStack.size()-1; i >= 0; --i)
   {
      Psi.push_back(RightStack[i]);
   }

   MatrixOperator M = CenterOp;
   ComponentType C = prod(LeftTop, M);
   M = TruncateBasis1(C);
   Psi.push_front(C);
   for (int i = LeftStack.size()-1; i >= 0; --i)
   {
      C = prod(LeftStack[i], M);
      M = TruncateBasis1(C);
      Psi.push_front(C);
   }

   // incorporate M back into the wavefunction
   Psi.set_front(prod(M, Psi.get_front()));
   return Psi;
}

struct WavefunctionCompat
{
   LinearWavefunction Psi;
};

inline
PStream::opstream& operator<<(PStream::opstream& out, WavefunctionCompat const& x)
{
   PANIC("Attempt to write a WavefunctionCompat");
   return out;
}

inline 
PStream::ipstream& operator>>(PStream::ipstream& in, WavefunctionCompat& x)
{
   x.Psi = LoadWavefunction_0_7_3(in);
   return in;
}

#endif
