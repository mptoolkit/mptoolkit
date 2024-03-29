// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/reducible.cc
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

namespace Tensor
{

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>::ReducibleTensor()
{
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>::ReducibleTensor(basis1_type const& Basis)
   : Basis1_(Basis), Basis2_(Basis)
{
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>::ReducibleTensor(basis1_type const& Basis1,
                                               basis2_type const& Basis2)
   : Basis1_(Basis1), Basis2_(Basis2)
{
}

template <typename T, typename B1, typename B2, typename S>
template <typename U, typename US>
ReducibleTensor<T, B1, B2, S>::ReducibleTensor(IrredTensor<U, B1, B2, US> const& x)
   : Basis1_(x.Basis1()), Basis2_(x.Basis2())
{
   data_[x.TransformsAs()] = x;
}

template <typename T, typename B1, typename B2, typename S>
template <typename U, typename US>
ReducibleTensor<T, B1, B2, S>::ReducibleTensor(ReducibleTensor<U, basis1_type, basis2_type, US> const& x)
   : Basis1_(x.Basis1()), Basis2_(x.Basis2()), data_(x.data_.begin(), x.data_.end())
{
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>&
ReducibleTensor<T, B1, B2, S>::operator+=(ReducibleTensor const& x)
{
   if (this->is_null())
   {
      (*this) = x;
      return *this;
   }
   for (typename ReducibleTensor<T,B1,B2,S>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      this->project(I->TransformsAs()) += *I;
   }
   return *this;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>&
ReducibleTensor<T, B1, B2, S>::operator-=(ReducibleTensor const& x)
{
   if (this->is_null())
   {
      (*this) = -x;
      return *this;
   }
   for (typename ReducibleTensor<T,B1,B2,S>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      this->project(I->TransformsAs()) -= *I;
   }
   return *this;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>&
ReducibleTensor<T, B1, B2, S>::operator+=(IrredTensorType const& x)
{
   if (this->is_null())
      (*this) = x;
   else
      this->project(x.TransformsAs()) += x;
   return *this;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>&
ReducibleTensor<T, B1, B2, S>::operator-=(IrredTensorType const& x)
{
   if (this->is_null())
      (*this) = -x;
   else
      this->project(x.TransformsAs()) -= x;
   return *this;
}

template <typename T, typename B1, typename B2, typename S>
template <typename U>
ReducibleTensor<T, B1, B2, S>&
ReducibleTensor<T, B1, B2, S>::operator*=(U const& x)
{
   for (iterator I = this->begin(); I != this->end(); ++I)
   {
      (*I) *= x;
   }
   return *this;
}

template <typename T>
struct TransformsAsEqualTo
{
   TransformsAsEqualTo(QuantumNumber const& q) : q_(q) {}

   typedef T argument_type;
   typedef bool result_type;

   bool operator()(T const& x) const
   {
      return x.TransformsAs() == q_;
   };

   QuantumNumber q_;
};

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T,B1,B2,S>
ReducibleTensor<T, B1, B2, S>::project(QuantumNumber const& q) const
{
   const_iterator I = std::find_if(this->begin(), this->end(),
                                   TransformsAsEqualTo<IrredTensor<T,B1,B2,S> >(q));
   if (I == this->end())
   {
      return IrredTensor<T,B1,B2,S>(Basis1_, Basis2_, q);
   }
   return *I;
}

template <typename T, typename B1, typename B2, typename S>
std::set<QuantumNumber>
ReducibleTensor<T, B1, B2, S>::components() const
{
   std::set<QuantumNumber> Result;
   for (const_map_iterator I = data_.begin(); I != data_.end(); ++I)
   {
      DEBUG_CHECK_EQUAL(I->first, I->second.TransformsAs())("Structure error in ReducibleTensor");
      Result.insert(I->first);
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T,B1,B2,S>&
ReducibleTensor<T, B1, B2, S>::project(QuantumNumber const& q)
{

   DEBUG_CHECK_EQUAL(q.GetSymmetryList(), this->GetSymmetryList());
   map_iterator I = data_.find(q);
   if (I == data_.end())
   {
      I = data_.insert(std::make_pair(q, IrredTensorType(Basis1_, Basis2_, q))).first;
   }
   return I->second;
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T,B1,B2,S> const&
ReducibleTensor<T, B1, B2, S>::project_assert(QuantumNumber const& q) const
{
   const_map_iterator I = data_.find(q);
   DEBUG_CHECK(I != data_.end());
   return I->second;
}

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T,B1,B2,S>
ReducibleTensor<T, B1, B2, S>::scalar() const
{
   return this->project(QuantumNumber(this->GetSymmetryList()));
}

template <typename T, typename B1, typename B2, typename S>
inline
IrredTensor<T,B1,B2,S>&
ReducibleTensor<T, B1, B2, S>::scalar()
{
   return this->project(QuantumNumber(this->GetSymmetryList()));
}

template <typename T, typename B1, typename B2, typename S>
void
ReducibleTensor<T, B1, B2, S>::check_structure() const
{
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      I->check_structure();
   }
}

template <typename T, typename B1, typename B2, typename S>
void
ReducibleTensor<T, B1, B2, S>::debug_check_structure() const
{
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      I->debug_check_structure();
   }
}

template <typename T, typename B1, typename B2, typename S>
PStream::opstream& operator<<(PStream::opstream& out, ReducibleTensor<T, B1, B2, S> const& Op)
{
   return out << Op.Basis1_ << Op.Basis2_ << Op.data_;
}

template <typename T, typename B1, typename B2, typename S>
PStream::ipstream& operator>>(PStream::ipstream& in, ReducibleTensor<T, B1, B2, S>& Op)
{
   return in >> Op.Basis1_ >> Op.Basis2_ >> Op.data_;
}

template <typename T, typename B1, typename B2, typename S>
std::ostream& operator<<(std::ostream& out, ReducibleTensor<T, B1, B2, S> const& Op)
{
   for (typename ReducibleTensor<T, B1, B2, S>::const_iterator I = Op.begin();
        I != Op.end(); ++I)
   {
      out << (*I) << '\n';
   }
   return out;
}

template <typename T, typename B1, typename B2, typename S>
std::string show_projections(ReducibleTensor<T, B1, B2, S> const& Op)
{
   std::string Result;
   for (typename ReducibleTensor<T, B1, B2, S>::const_iterator I = Op.begin();
        I != Op.end(); ++I)
   {
      Result += "Components that transform as "
         + boost::lexical_cast<std::string>(I->first)
         + ":\n";
      Result += show_projections(I->second);
   }
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
ReducibleTensor<T, B1, B2, S>
CoerceSymmetryList(ReducibleTensor<T, B1, B2, S> const& t, SymmetryList const& sl)
{
   ReducibleTensor<T, B1, B2, S> Result(t);
   CoerceSymmetryListInPlace(Result, sl);
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
void
CoerceSymmetryListInPlace(ReducibleTensor<T, B1, B2, S>& t, SymmetryList const& sl)
{
   CoerceSymmetryListInplace(t.Basis1_, sl);
   CoerceSymmetryListInplace(t.Basis2_, sl);
   for (auto& x : t.data_)
   {
      CoerceSymmetryListInplace(x.first, sl);
      CoerceSymmetryListInplace(x.second, sl);
   }
}

} // namespace Tensor
