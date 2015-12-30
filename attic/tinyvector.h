// -*- C++ -*- $Id$

/*
// -*- mode:C++; tab-width:4; c-basic-offset:4; indent-tabs-mode:nil -*-

  tinyvector.h

  Ian McCulloch
  created 29/10/98

  Renamed TinyVector, and modified 2000-06-04.  No longer uses exprtemplate.h,
  TinyVecImp removed.

  The TinyVec class implements a vector of objects, where the number of objects is
  a (small) compile-time constant.  TinyVecImp is used as a base class, that defines
  various functions.  The derived class TinyVec<class T, int N> is specialized
  for small values of N, to provide a constructor to initialize all elements
  separately.  For larger values of N (currently > 4), the default TinyVec class applies, 
  and there is no such constructor.  TinyVecImp should not be used outside this header.

  TinyVec<T,0> is a legitimate class.  It contains no members.

*/

#if !defined(TINYVECTOR_H_G675F756H4w5J78975H567)
#define TINYVECTOR_H_G675F756H4w5J78975H567

#include <iostream>
#include <iomanip>
#include "trace.h"
#include "config.h"
#include "ctassert.h"

namespace TinyAlgebra
{

//using namespace PStream;

template <class T, int N> class TinyVector;

template <class T, int N> 
class TinyVector
{
   public:
      typedef T                 value_type;
      typedef T* restrict       pointer;
      typedef T const* restrict const_pointer;
      typedef T& restrict       reference;
      typedef T const& restrict const_reference;
      typedef T* restrict       iterator;
      typedef T const* restrict const_iterator;

      static int const size = N;

      TinyVector() {}
      ~TinyVector() {}

      TinyVector(const TinyVector<T, N>& Vec);

      TinyVector<T,N>& operator=(const TinyVector<T,N>& Vec);

      iterator begin()             { return Data; }
      const_iterator begin() const { return Data; }

      iterator end()             { return Data + size; }
      const_iterator end() const { return Data + size; }

      pointer data()             { return Data; }
      const_pointer data() const { return Data; }

      reference operator[](int i);
      value_type operator[](int i) const;

      // 1 arg initial value ctor is explicit.  TinyVector<T, 1> is partially specialized
      explicit TinyVector(T InitValue);

      TinyVector(T x0, T x1) 
      {
         typedef ct_assert<size == 2> size_check;
         Data[0] = x0;
         Data[1] = x1;
      }

      TinyVector(T x0, T x1, T x2) 
      {
         typedef ct_assert<size == 3> size_check;
         Data[0] = x0;
         Data[1] = x1;
         Data[2] = x2;
      }

      TinyVector(T x0, T x1, T x2, T x3) 
      {
         typedef ct_assert<size == 4> size_check;
         Data[0] = x0;
         Data[1] = x1;
         Data[2] = x2;
         Data[3] = x3;
      }
      // need to add more initializers as they are needed

      TinyVector<T, N>& operator*=(T t)
      {
         for (int i = 0; i < N; ++i)
	 {
            Data[i] *= t;
	 }
         return *this;
      }

      TinyVector<T, N>& operator/=(T t)
      {
         for (int i = 0; i < N; ++i)
	 {
            Data[i] /= t;
	 }
         return *this;
      }

      TinyVector<T, N>& operator+=(TinyVector<T, N> const& v)
      {
         for (int i = 0; i < N; ++i)
	 {
            Data[i] += v.Data[i];
	 }
         return *this;
      }

      TinyVector<T, N>& operator-=(TinyVector<T, N> const& v)
      {
         for (int i = 0; i < N; ++i)
	 {
            Data[i] -= v.Data[i];
	 }
         return *this;
      }

   private:
      bool in_range(int i) const { return i >= 0 && i < size; }

      value_type Data[size];      
};

template <class T>
class TinyVector<T, 1>
{
   public:
      typedef T                 value_type;
      typedef T* restrict       pointer;
      typedef T const* restrict const_pointer;
      typedef T& restrict       reference;
      typedef T const& restrict const_reference;
      typedef T* restrict       iterator;
      typedef T const* restrict const_iterator;

      static int const size = 1;

      TinyVector() { }
      TinyVector(TinyVector<T, size> const& Vec) { Data = Vec.Data; }
      TinyVector(T x0) { Data = x0; }

      TinyVector<T,1>& operator=(TinyVector<T, size> const& Vec) { Data = Vec.Data; return *this; }

      iterator begin()             { return Data; }
      const_iterator begin() const { return Data; }

      iterator end()             { return Data + size; }
      const_iterator end() const { return Data + size; }

      pointer data()             { return Data; }
      const_pointer data() const { return Data; }

      reference operator[](int i);
      value_type operator[](int i) const;

      TinyVector<T, 1>& operator*=(T t)
      {
	 Data *= t;
	 return *this;
      }

      TinyVector<T, 1>& operator/=(T t)
      {
	 Data /= t;
         return *this;
      }

      TinyVector<T, 1>& operator+=(TinyVector<T, 1> const& v)
      {
	 Data += v.Data;
	 return *this;
      }

      TinyVector<T, 1>& operator-=(TinyVector<T, 1> const& v)
      {
	 Data -= v.Data;
	 return *this;
      }

   private:
      bool in_range(int i) const { return i == 0; }

      value_type Data;      
};

// For completeness, we have TinyVector<T, 0>

template <class T>
class TinyVector<T, 0>
{
   public:
      typedef T                 value_type;
      typedef T* restrict       pointer;
      typedef T const* restrict const_pointer;
      typedef T& restrict       reference;
      typedef T const& restrict const_reference;
      typedef T* restrict       iterator;
      typedef T const* restrict const_iterator;

      static int const size = 0;

      TinyVector() {}
      TinyVector(TinyVector<T, size> const& Vec) {}
      // we still have an initial value ctor
      explicit TinyVector(T InitValue) {}

      TinyVector<T,0>& operator=(TinyVector<T, size> const& Vec) { return *this; }

      iterator begin()             { return NULL; }
      const_iterator begin() const { return NULL; }

      iterator end()             { return NULL; }
      const_iterator end() const { return NULL; }

      pointer data()             { return NULL; }
      const_pointer data() const { return NULL; }

      // no operator[], not meaningful for 0 length vector

      // we probably should have arithmetic operators though
};

template <class T>
inline
bool operator<(TinyVector<T, 1> const& x, TinyVector<T, 1> const& y)
{
   return x[0] < y[0];
}

template <class T>
inline
bool operator<(TinyVector<T, 2> const& x, TinyVector<T, 2> const& y)
{
   return x[0] < y[0] || (x[0] == y[0] && x[1] < y[1]);
}

template <class T, int N>
std::ostream& operator<<(std::ostream& out, TinyVector<T, N> const& Vec)
{
   out << '(';
   for (int i = 0; i < N; ++i)
   {
      if (i != 0) out << ',';
      out << Vec[i];
   }
   out << ')';
   return out;
}      

template <class T, int N>
std::istream& operator>>(std::istream& in, TinyVector<T, N>& Vec)
{
   char junk;
   in >> std::skipws >> junk;
   CHECK(junk == '(');
   //   in >> junk;
   for (int i = 0; i < N; ++i)
   {
      //      if (i != 0) in >> std::skipws >> junk;
      if (i != 0) in >> junk;
      in >> Vec[i];
   }
   //   in >> std::skipws >> junk;
   in >> junk;
   CHECK(junk == ')');
   return in;
}      

template <class T, int N>
T dot(TinyVector<T, N> const& x, TinyVector<T, N> const& y)
{
   T Result = 0;
   for (int i = 0; i < N; ++i)
   {
      Result += x[i] * y[i];
   }
   return Result;
}

template <class T, int N>
T norm(TinyVector<T, N> const& x)
{
   T Result = 0;
   for (int i = 0; i < N; ++i)
   {
      Result += x[i] * x[i];
   }
   return Result;
}

// inlines

template <class T, int N>
inline
TinyVector<T, N>::TinyVector(const TinyVector<T, N>& Vec)
{
   for (int i = 0; i < size; ++i)
      Data[i] = Vec.Data[i];
}

template <class T, int N>
inline
TinyVector<T,N>& 
TinyVector<T, N>::operator=(const TinyVector<T,N>& Vec)
{
   for (int i = 0; i < size; ++i)
      Data[i] = Vec.Data[i];

   return *this;
}

template <class T, int N>
inline
TinyVector<T, N>::TinyVector(T InitValue)
{
   for (int i = 0; i < size; ++i)
      Data[i] = InitValue;
}

template <class T, int N>
inline
TinyVector<T, N>::reference 
TinyVector<T, N>::operator[](int i)
{
   PRECONDITION(in_range(i));
   return Data[i];
}

template <class T, int N>
inline
T
TinyVector<T, N>::operator[](int i) const
{
   PRECONDITION(in_range(i));
   return Data[i];
}

template <class T>
inline
TinyVector<T, 1>::reference 
TinyVector<T, 1>::operator[](int i)
{
   PRECONDITION(in_range(i));
   return Data[i];
}

template <class T>
inline
T
TinyVector<T, 1>::operator[](int i) const
{
   PRECONDITION(in_range(i));
   return Data[i];
}

} // namespace TinyAlgebra

#endif // sentry

