// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/centerwavefunction.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "centerwavefunction.h"

CenterWavefunction::CenterWavefunction(LinearWavefunction const& Psi)
   : Attr(Psi.Attributes())
{
   LinearWavefunction::const_iterator I = Psi.begin();
   LeftTop = *I;
   CenterOp = ExpandBasis2(LeftTop);
   I = Psi.end();
   --I;
   while (I != Psi.begin())
   {
      RightStack.push_back(new ComponentType(*I));
      --I;
   }
   RightTop = *RightStack.back().load();
   RightStack.pop_back();
}

CenterWavefunction::CenterWavefunction(MatrixOperator const& C, LinearWavefunction const& Psi)
   : Attr(Psi.Attributes())
{
   LinearWavefunction::const_iterator I = Psi.begin();
   LeftTop = prod(C, *I);
   CenterOp = ExpandBasis2(LeftTop);
   I = Psi.end();
   --I;
   while (I != Psi.begin())
   {
      RightStack.push_back(new ComponentType(*I));
      --I;
   }
   RightTop = *RightStack.back().load();
   RightStack.pop_back();
}

VectorBasis CenterWavefunction::LeftVacuumBasis() const
{
   return this->LeftSize() == 0 ? this->Center().Basis1() : this->LookupLeft(0).Basis1();
}

VectorBasis CenterWavefunction::RightVacuumBasis() const
{
   return this->RightSize() == 0 ? this->Center().Basis2() : this->LookupRight(0).Basis2();
}

CenterWavefunction& CenterWavefunction::normalize()
{
   *this *= 1.0 / norm_2(*this); return *this;
}

CenterWavefunction& CenterWavefunction::Normalize()
{
   *this *= 1.0 / norm_2(*this); return *this;
}

QuantumNumber CenterWavefunction::TransformsAs() const
{
   // a null matrix product state is assumed to transform as a scalar
   if (this->LookupLeft(0).Basis1().size() == 0)
      return QuantumNumber(this->Center().GetSymmetryList());

   CHECK_EQUAL(this->LookupLeft(0).Basis1().size(), 1);
   return this->LookupLeft(0).Basis1()[0];
}

void CenterWavefunction::RotateLeft()
{
   // DEBUG_PRECONDITION(this->Right() satisfies orthogonality constraint)

   // Shift the Center matrix to the left to become the Pivot,
   // and apply the orthogonality constraint (singular value decomposition)
   this->PushRight(prod(this->Left(), this->Center()));
   this->PopLeft();
   this->Center() = TruncateBasis1(this->Right());
}

void CenterWavefunction::RotateRight()
{
   // DEBUG_PRECONDITION(this->Left() satisfies orthogonality constraint)

   // Shift the Center matrix to the right to become the Pivot,
   // and apply the orthogonality constraint (singular value decomposition)
   this->PushLeft(prod(this->Center(), this->Right()));
   this->PopRight();
   this->Center() = TruncateBasis2(this->Left());
}

void CenterWavefunction::RotateLeftExpand()
{
   this->PushRight(prod(this->Left(), this->Center()));
   this->PopLeft();
   this->Center() = ExpandBasis1(this->Right());
}

void CenterWavefunction::RotateRightExpand()
{
   this->PushLeft(prod(this->Center(), this->Right()));
   this->PopRight();
   this->Center() = ExpandBasis2(this->Left());
}

#if 0
void CenterWavefunction::RotateLeftTruncate(int MaxStates)
{
   // DEBUG_PRECONDITION(this->Right() satisfies orthogonality constraint)

   // Shift the Center matrix to the left to become the Pivot,
   // and apply the orthogonality constraint (singular value decomposition)
   this->PushRight(prod(this->Left(), this->Center()));
   this->PopLeft();
   this->Center() = TruncateBasis1(this->Right(), MaxStates);
}

void CenterWavefunction::RotateRightTruncate(int MaxStates)
{
   // DEBUG_PRECONDITION(this->Left() satisfies orthogonality constraint)

   // Shift the Center matrix to the right to become the Pivot,
   // and apply the orthogonality constraint (singular value decomposition)
   this->PushLeft(prod(this->Center(), this->Right()));
   this->PopRight();
   this->Center() = TruncateBasis2(this->Left(), MaxStates);
}
#endif

LinearWavefunction
CenterWavefunction::AsLinearWavefunction() const
{
   LinearWavefunction Psi(this->GetSymmetryList());
   Psi.Attributes() = this->Attributes();
   for (int i = this->RightSize()-1; i >= 0; --i)
   {
      Psi.push_back(this->LookupRight(i));
   }

   MatrixOperator M = this->Center();
   for (int i = this->LeftSize()-1; i >= 0; --i)
   {
      ComponentType C = prod(this->LookupLeft(i), M);
      M = TruncateBasis1(C);
      Psi.push_front(C);
   }

   // incorporate M back into the wavefunction
   Psi.set_front(prod(M, Psi.get_front()));
   return Psi;
}

void CenterWavefunction::RotateToNormalForm()
{
   while (this->LeftSize() > 1)
      this->RotateLeft();
}

CenterWavefunction::CenterWavefunction()
{
}

SymmetryList CenterWavefunction::GetSymmetryList() const
{
   if (this->Left().is_null())
      return this->Right().GetSymmetryList();
   else
      return this->Left().GetSymmetryList();
}


int CenterWavefunction::LeftSize() const
{
   if (LeftTop.is_null()) return 0;
   return LeftStack.size()+1;
}


int CenterWavefunction::RightSize() const
{
   if (RightTop.is_null()) return 0;
   return RightStack.size()+1;
}


void CenterWavefunction::PushLeft(MPStateComponent const& L)
{
   PRECONDITION(LeftTop.is_null() || LeftTop.GetSymmetryList() == L.GetSymmetryList());
   PRECONDITION(RightTop.is_null() || RightTop.GetSymmetryList() == L.GetSymmetryList());
   if (!LeftTop.is_null()) LeftStack.push_back(new MPStateComponent(LeftTop));
   LeftTop = L;
}


void CenterWavefunction::PushRight(MPStateComponent const& R)
{
   PRECONDITION(LeftTop.is_null() || LeftTop.GetSymmetryList() == R.GetSymmetryList());
   PRECONDITION(RightTop.is_null() || RightTop.GetSymmetryList() == R.GetSymmetryList());
   if (!RightTop.is_null()) RightStack.push_back(new MPStateComponent(RightTop));
   RightTop = R;
}


void CenterWavefunction::PopLeft()
{
   PRECONDITION(!LeftTop.is_null());
   if (LeftStack.empty())
   {
      LeftTop = ComponentType();
   }
   else
   {
      LeftTop = *LeftStack.back().load();
      LeftStack.pop_back();
   }
}


void CenterWavefunction::PopRight()
{
   PRECONDITION(!RightTop.is_null());
   if (RightStack.empty())
   {
      RightTop = ComponentType();
   }
   else
   {
      RightTop = *RightStack.back().load();
      RightStack.pop_back();
   }
}
