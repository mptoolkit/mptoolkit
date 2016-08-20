// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/mapvector.cc
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

namespace LinearAlgebra
{

template <typename T>
template <typename U>
inline
MapVector<T>::MapVector(U const& x, typename boost::enable_if<is_vector<U> >::type*)
   : Size_(Size<U>()(x))
{
   assign(*this, x);
}

template <typename T>
template <typename U>
inline
MapVector<T>::MapVector(NoAliasProxy<U> const& x, typename boost::enable_if<is_vector<U> >::type*)
   : Size_(Size<U>()(x))
{
   assign(*this, x.value());
}

template <typename T>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, MapVector<T>&>::type
MapVector<T>::operator=(U const& x)
{
   MapVector<T> Temp(x);
   this->swap(Temp);
   return *this;
}

template <typename T>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, MapVector<T>&>::type
MapVector<T>::operator=(NoAliasProxy<U> const& x)
{
   Data_.clear();
   Size_ = Size<U>(x.value());
   assign(*this, x.value());
   return *this;
}

template <typename T>
inline
T const&
MapVector<T>::operator[](size_type n) const
{
   base_const_iterator I = Data_.find(n);
   return (I == Data_.end()) ? static_zero_or_die<T>() : I->second;
}

template <typename T>
inline
T&
MapVector<T>::operator[](size_type n)
{
   base_iterator I = Data_.find(n);
   if (I == Data_.end()) I = Data_.insert(std::pair<size_type, T>(n, zero_or_default<T>())).first;
   return I->second;
}

template <typename T>
template <typename U>
inline
void
MapVector<T>::set_element(size_type n, U const& x)
{
   base_iterator I = Data_.find(n);
   if (I != Data_.end())
      I->second = x;
   else
      Data_.insert(std::pair<size_type, T>(n, x));
}

template <typename T>
template <typename U>
inline
void
MapVector<T>::add_element(size_type n, U const& x)
{
   base_iterator I = Data_.find(n);
   if (I != Data_.end())
      I->second += x;
   else
      Data_.insert(std::pair<size_type, T>(n, x));
}

template <typename T>
template <typename U, typename Float>
inline
void
MapVector<T>::add_element_cull(size_type n, U const& x, Float const& Tol)
{
   // Early return if the element is zero, eg it might be a value_with_zero<T>
   if (LinearAlgebra::is_zero(x))
      return;

   T Temp(x);
   base_iterator I = Data_.find(n);
   if (I != Data_.end())
   {
      auto NormSq = norm_err_sq(Temp) + norm_err_sq(I->second);
      I->second += Temp;
      if (norm_err_sq(I->second) < Tol*Tol * NormSq)
      {
         Data_.erase(I);
      }
   }
   else
      Data_.insert(std::pair<size_type, T>(n, Temp));
}


template <typename T>
template <typename U>
inline
void
MapVector<T>::subtract_element(size_type n, U const& x)
{
   base_iterator I = Data_.find(n);
   if (I != Data_.end())
      I->second -= x;
   else
      Data_.insert(std::pair<size_type, T>(n, -x));
}

template <typename T>
template <typename U, typename Float>
inline
void
MapVector<T>::subtract_element_cull(size_type n, U const& x, Float const& Tol)
{
   if (LinearAlgebra::is_zero(x))
      return;

   T Temp(x);
   base_iterator I = Data_.find(n);
   if (I != Data_.end())
   {
      auto NormSq = norm_err_sq(Temp) + norm_err_sq(I->second);
      I->second -= Temp;
      if (norm_err_sq(I->second) < Tol*Tol * NormSq)
      {
         Data_.erase(I);
      }
   }
   else
      Data_.insert(std::pair<size_type, T>(n, -Temp));
}

template <typename T>
inline
void
MapVector<T>::zero_element(size_type n)
{
   Data_.erase(n);
}

} // namespace LinearAlgebra
