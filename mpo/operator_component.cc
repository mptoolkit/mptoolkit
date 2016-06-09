// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_component.cc
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

inline
OperatorComponent& 
OperatorComponent::operator+=(OperatorComponent const& x)
{
   DEBUG_CHECK_EQUAL(LocalBasis1_, x.LocalBasis1_);
   DEBUG_CHECK_EQUAL(LocalBasis2_, x.LocalBasis2_);
   DEBUG_CHECK_EQUAL(Basis1_, x.Basis1_);
   DEBUG_CHECK_EQUAL(Basis2_, x.Basis2_);
   Data_ += x.Data_;
   return *this;
}

inline
OperatorComponent& 
OperatorComponent::operator-=(OperatorComponent const& x)
{
   DEBUG_CHECK_EQUAL(LocalBasis1_, x.LocalBasis1_);
   DEBUG_CHECK_EQUAL(LocalBasis2_, x.LocalBasis2_);
   DEBUG_CHECK_EQUAL(Basis1_, x.Basis1_);
   DEBUG_CHECK_EQUAL(Basis2_, x.Basis2_);
   Data_ -= x.Data_;
   return *this;
}

inline
OperatorComponent::value_type
OperatorComponent::operator()(int i, int j) const
{
   const_inner_iterator I = LinearAlgebra::iterate_at(Data_, i,j);
   if (I)
      return *I;

   else return value_type(LocalBasis1_, LocalBasis2_);
}

inline
OperatorComponent::value_type&
OperatorComponent::operator()(int i, int j)
{
   inner_iterator I = LinearAlgebra::iterate_at(Data_, i,j);
   if (!I)
   {
      Data_(i,j) = value_type(LocalBasis1_, LocalBasis2_);
      I = LinearAlgebra::iterate_at(Data_, i,j);
      DEBUG_CHECK(!!I);
   }
   return *I;
}

inline
OperatorComponent::value_type
OperatorComponent::top_left() const
{
   return this->operator()(0,0);
}

inline
OperatorComponent::value_type
OperatorComponent::bottom_right() const
{
   return this->operator()(this->size1()-1, this->size2()-1);
}

inline
bool
OperatorComponent::is_lower_triangular() const
{
   for (const_iterator I = iterate(Data_); I; ++I)
   {
      for (const_inner_iterator J = iterate(I); J; ++J)
      {
	 if (J.index1() > J.index2())
	    return false;
      }
   }
   return true;
}

inline
void
OperatorComponent::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}
