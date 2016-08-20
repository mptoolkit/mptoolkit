// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/splitoperator.h
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
// SplitOperator - backwards compatibility wrapper for
// MPOperator.
//

#if !defined(MPOPERATOR_COMPAT_H_ASDKJFHSDKJ23489753829)
#define MPOPERATOR_COMPAT_H_ASDKJFHSDKJ23489753829

#include "mpoperator.h"
#include <algorithm>

class SplitOperator
{
   public:
      typedef MPOpComponent ComponentType;
      typedef SimpleOperator OperatorType;

      SplitOperator() {}

      SplitOperator(MPOperator const& Op);

      SplitOperator(SplitOperator const& Op);
      SplitOperator& operator=(SplitOperator const& Op);

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Op_.GetSymmetryList(); }

      QuantumNumbers::QuantumNumber TransformsAs() const { return Op_.TransformsAs(); }

      int size() const;
      int LeftSize() const;
      int RightSize() const;

      ComponentType Left() const;
      ComponentType Right() const;
      OperatorType Center() const;

      ComponentType LookupLeft(int n) const;
      ComponentType LookupRight(int n) const;

      BasisList LeftVacuumBasis() const;
      BasisList RightVacuumBasis() const;

      SplitOperator& operator*=(double x);
      SplitOperator& operator*=(std::complex<double> const& x);

      SplitOperator& operator+=(SplitOperator const& Op);
      SplitOperator& operator-=(SplitOperator const& Op);

      void RotateLeft();
      void RotateRight();
      void RotateToNormalForm();

      MPOperator const& AsMPOperator() const { return Op_; }

   private:
      void ResetIterFromLoc();

      MPOperator Op_;
      int Loc_;
      MPOperator::const_iterator I_;

   friend PStream::opstream& operator<<(PStream::opstream& out, SplitOperator const& Op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SplitOperator& Op);
};

inline
PStream::opstream& operator<<(PStream::opstream& out, SplitOperator const& Op)
{
   out << Op.Op_;
   out << Op.Loc_;
   return out;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, SplitOperator& Op)
{
   in >> Op.Op_;
   in >> Op.Loc_;
   Op.ResetIterFromLoc();
   return in;
}

inline
SplitOperator::SplitOperator(MPOperator const& Op)
   : Op_(Op), Loc_(1), I_(Op_.begin())
{
}

inline
SplitOperator::SplitOperator(SplitOperator const& Op)
   : Op_(Op.Op_), Loc_(Op.Loc_)
{
   this->ResetIterFromLoc();
}

inline
SplitOperator& SplitOperator::operator=(SplitOperator const& Op)
{
   Op_ = Op.Op_;
   Loc_ = Op.Loc_;
   this->ResetIterFromLoc();
   return *this;
}

inline
int SplitOperator::size() const
{
   return Op_.size();
}

inline
int SplitOperator::LeftSize() const
{
   return Loc_;
}

inline
int SplitOperator::RightSize() const
{
   return Op_.size() - Loc_;
}

inline
SplitOperator::ComponentType
SplitOperator::Left() const
{
   return *I_;
}

inline
SplitOperator::ComponentType
SplitOperator::Right() const
{
   MPOperator::const_iterator J = I_;
   ++J;
   return *J;
}

inline
SplitOperator::OperatorType
SplitOperator::Center() const
{
   return OperatorType::make_identity(I_->Basis2());
}

inline
SplitOperator::ComponentType SplitOperator::LookupLeft(int n) const
{
   MPOperator::const_iterator J = Op_.begin();
   std::advance(J, n);
   return *J;
}

inline
SplitOperator::ComponentType SplitOperator::LookupRight(int n) const
{
   MPOperator::const_iterator J = Op_.end();
   std::advance(J, -n);
   --J;
   return *J;
}

inline
BasisList SplitOperator::LeftVacuumBasis() const
{
   return Op_.Basis1();
}

inline
BasisList SplitOperator::RightVacuumBasis() const
{
   return Op_.Basis2();
}

inline
SplitOperator& SplitOperator::operator*=(double x)
{
   Op_ *= x;
   return *this;
}

inline
SplitOperator&
SplitOperator::operator*=(std::complex<double> const& x)
{
   Op_ *= x;
   return *this;
}

inline
SplitOperator& SplitOperator::operator+=(SplitOperator const& Op)
{
   Op_ += Op.AsMPOperator();
   this->ResetIterFromLoc();
   return *this;
}

inline
SplitOperator& SplitOperator::operator-=(SplitOperator const& Op)
{
   Op_ -= Op.AsMPOperator();
   this->ResetIterFromLoc();
   return *this;
}

inline
void SplitOperator::RotateLeft()
{
   CHECK(Loc_ > 1);
   --Loc_;
   --I_;
}

inline
void SplitOperator::RotateRight()
{
   CHECK(Loc_ < Op_.size());
   ++Loc_;
   ++I_;
}

inline
void SplitOperator::RotateToNormalForm()
{
   Loc_ = 1;
   I_ = Op_.begin();
}

inline
void SplitOperator::ResetIterFromLoc()
{
   I_ = Op_.begin();
   std::advance(I_, Loc_-1);
}

inline
SplitOperator prod(SplitOperator const& x, SplitOperator const& y, QuantumNumbers::QuantumNumber const& q)
{
   return SplitOperator(prod(x.AsMPOperator(), y.AsMPOperator(), q));
}

inline
SplitOperator prod(SplitOperator const& x, SplitOperator const& y)
{
   return SplitOperator(prod(x.AsMPOperator(), y.AsMPOperator()));
}

#endif
