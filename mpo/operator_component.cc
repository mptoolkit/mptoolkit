// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator_component.cc
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

template <typename T>
BasicOperatorComponent<T>::BasicOperatorComponent(BasisList const& LocalB,
						  BasisList const& B1, BasisList const& B2)
   : LocalBasis1_(LocalB), LocalBasis2_(LocalB), Basis1_(B1), Basis2_(B2),
     Data_(Basis1_.size(), Basis2_.size())
{
}

template <typename T>
BasicOperatorComponent<T>::BasicOperatorComponent(BasisList const& LocalB1, BasisList const& LocalB2,
						  BasisList const& B1, BasisList const& B2)
   : LocalBasis1_(LocalB1), LocalBasis2_(LocalB2), Basis1_(B1), Basis2_(B2),
     Data_(Basis1_.size(), Basis2_.size())
{
}


template <typename T>
void
BasicOperatorComponent<T>::check_structure() const
{
   for (auto const& r : (*this))
   {
      for (auto const& c : r)
      {
	 for (auto const& q : c.value)
	 {
	    CHECK(is_transform_target(this->qn2(c.col()), q.TransformsAs(), this->qn1(r.row())))
	       (this->qn2(c.col()))(q.TransformsAs())(this->qn1(r.row()));
	 }
      }
   }
}

template <typename T>
inline
BasicOperatorComponent<T>&
BasicOperatorComponent<T>::operator+=(BasicOperatorComponent<T> const& x)
{
   DEBUG_CHECK_EQUAL(LocalBasis1_, x.LocalBasis1_);
   DEBUG_CHECK_EQUAL(LocalBasis2_, x.LocalBasis2_);
   DEBUG_CHECK_EQUAL(Basis1_, x.Basis1_);
   DEBUG_CHECK_EQUAL(Basis2_, x.Basis2_);
   Data_ += x.Data_;
   return *this;
}

template <typename T>
inline
BasicOperatorComponent<T>&
BasicOperatorComponent<T>::operator-=(BasicOperatorComponent<T> const& x)
{
   DEBUG_CHECK_EQUAL(LocalBasis1_, x.LocalBasis1_);
   DEBUG_CHECK_EQUAL(LocalBasis2_, x.LocalBasis2_);
   DEBUG_CHECK_EQUAL(Basis1_, x.Basis1_);
   DEBUG_CHECK_EQUAL(Basis2_, x.Basis2_);
   Data_ -= x.Data_;
   return *this;
}

#if 0
template <typename T>
inline
BasicOperatorComponent<T>::value_type
BasicOperatorComponent<T>::operator()(int i, int j) const
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
#endif

template <typename T>
inline
typename BasicOperatorComponent<T>::value_type
BasicOperatorComponent<T>::top_left() const
{
   typename data_type::row_type::const_iterator I = Data_[0].find(0);
   if (I == Data_[0].end())
      return value_type(LocalBasis1_, LocalBasis2_);
   // else
   return I.value;
}

template <typename T>
inline
typename BasicOperatorComponent<T>::value_type
BasicOperatorComponent<T>::bottom_right() const
{
   typename data_type::row_type::const_iterator I = Data_.back().find(Basis2_.size()-1);
   if (I == Data_.back().end())
      return value_type(LocalBasis1_, LocalBasis2_);
   // else
   return I.value;
}

template <typename T>
inline
void
BasicOperatorComponent<T>::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

template <typename T>
void print_structure(BasicOperatorComponent<T> const& Op, std::ostream& out, double UnityEpsilon)
{
   for (auto const& r : Op)
   {
      out << '[';
      int next_r = 0;
      for (auto const& c : r)
      {
	 while (next_r++ < r)
	    out << ' ';

	 if (c.value.size() > 1)
	    out << 'x';       // some compound operator
	 else if (!is_scalar(c.value))
	    out << 'v';       // a non-scalar
	 else
	 {
	    SimpleOperator Y = c.value.scalar();

	    complex x = PropIdent(Y);
	    if (x == 0.0)
	    {
	       complex x = PropIdent(scalar_prod(herm(Y),Y), UnityEpsilon);
	       if (norm_frob(x-1.0) < 1E-12)
		  out << 'U';      // a unitary
	       else
		  out << 's';      // a generic scalar
	    }
	    else if (norm_frob(x-1.0) < 1E-12)
	    {
	       out << 'I';         // the identity
	    }
	    else
	       out << 'i';         // something proportional to the identity
	 }
      }
      out << "]\n";
   }
}

template <typename T>
PStream::opstream&
operator<<(PStream::opstream& out, BasicOperatorComponent<T> const& Op)
{
   out << Op.LocalBasis1_ << Op.LocalBasis2_
       << Op.Basis1_ << Op.Basis2_
       << Op.Data_
      ;
   return out;
}

template <typename T>
PStream::ipstream&
operator>>(PStream::ipstream& in, BasicOperatorComponent<T>& Op)
{
   in >> Op.LocalBasis1_ >> Op.LocalBasis2_
      >> Op.Basis1_ >> Op.Basis2_
      >> Op.Data_
      ;
   return in;
}

template <typename T>
BasicOperatorComponent<T>
operator+(BasicOperatorComponent<T> const& A, BasicOperatorComponent<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());
   DEBUG_CHECK_EQUAL(A.Basis2(), B.Basis2());
   BasicOperatorComponent<T> Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), A.Basis2());
   Result.data() = A.data() + B.data();
   return Result;
}

template <typename T>
BasicOperatorComponent<T>
operator-(BasicOperatorComponent<T> const& A, BasicOperatorComponent<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.LocalBasis1(), B.LocalBasis1());
   DEBUG_CHECK_EQUAL(A.LocalBasis2(), B.LocalBasis2());
   DEBUG_CHECK_EQUAL(A.Basis1(), B.Basis1());
   DEBUG_CHECK_EQUAL(A.Basis2(), B.Basis2());
   BasicOperatorComponent<T> Result(A.LocalBasis1(), A.LocalBasis2(), A.Basis1(), A.Basis2());
   Result.data() = A.data() - B.data();
   return Result;
}
