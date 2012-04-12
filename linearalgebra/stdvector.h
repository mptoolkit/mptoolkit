/* -*- C++ -*- $Id$

  stdvector.h

  Interface declarations to enable expressions on std::vector.
 
  Created 2005-01-06 Ian McCulloch

  TODO: add allocator parameter
*/

#if !defined(STDVECTOR_H_DJKHCJKHTUIRYTORUL)
#define STDVECTOR_H_DJKHCJKHTUIRYTORUL

#include "vectorinterface.h"
#include "vecptriterator.h"
#include "vectormemproxy.h"
#include "slice.h"

namespace LinearAlgebra
{

// std::vector is a ContiguousVector

template <typename T>
struct interface<std::vector<T> >
{
   typedef CONTIGUOUS_VECTOR(T, std::vector<T>) type;
   typedef T value_type;
};

//
// iterators
//

template <typename T>
struct Iterate<std::vector<T> >
{
   typedef std::vector<T> const& argument_type;
   typedef VecPtrIterator<T const> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(&x[0], x.size(), 0);
   }
};

template <typename T>
struct Iterate<std::vector<T>&>
{
   typedef std::vector<T>& argument_type;
   typedef VecPtrIterator<T> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(&x[0], x.size(), 0);
   }
};

// resize

template <typename T>
struct Resize<std::vector<T>&>
{
   typedef void result_type;
   typedef std::vector<T>& first_argument_type;
   typedef size_type second_argument_type;
   void operator()(std::vector<T>& v, size_type n) const
   {
      v.resize(n);
   }
};

//
// Data
//

template <typename T>
struct Data<std::vector<T> >
{
   typedef T const* result_type;
   typedef std::vector<T> argument_type;
   result_type operator()(argument_type const& v) const { return &v[0]; }
};

template <typename T>
struct Data<std::vector<T>&>
{
   typedef T* result_type;
   typedef std::vector<T>& argument_type;
   result_type operator()(argument_type v) const { return &v[0]; }
};

// make_vec overload, useful for std::vector as an lvalue (eg, for expressions)

// TODO: this doesn't allow resizing of the left-hand-side.  Need a proper
// resizing proxy to allow this.
template <typename T>
inline
VectorMemProxy<T> make_vec(std::vector<T>& v)
{
   return VectorMemProxy<T>(&v[0], v.size());
}

template <typename T>
inline
VectorMemProxy<T const> make_vec(std::vector<T> const& v)
{
   return VectorMemProxy<T const>(&v[0], v.size());
}

// for convenience, we have an overload for noalias too

template <typename T>
inline
NoAliasProxy<VectorMemProxy<T> > noalias(std::vector<T>& v)
{
   return NoAliasProxy<VectorMemProxy<T> >(VectorMemProxy<T>(&v[0], v.size()));
}

template <typename T>
inline
NoAliasProxy<VectorMemProxy<T const> > noalias(std::vector<T> const& v)
{
   return NoAliasProxy<VectorMemProxy<T const> >(VectorMemProxy<T const>(&v[0], v.size()));
}

// operator+=, operator-=

template <typename Scalar, typename U>
inline
std::vector<Scalar>&
operator+=(std::vector<Scalar>& x, U const& y)
{
   // TODO: the choice of std::vector as temporary type isn't necessarily optimal here
   std::vector<Scalar> Temp(y);
   add(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
std::vector<Scalar>&
operator+=(std::vector<Scalar>& x, NoAliasProxy<U> const& y)
{
   add(x, y.value());
   return x;
}

template <typename Scalar, typename U>
inline
std::vector<Scalar>&
operator-=(std::vector<Scalar>& x, U const& y)
{
   // TODO: the choice of std::vector as temporary type isn't necessarily optimal here
   std::vector<Scalar> Temp(y);
   subtract(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
std::vector<Scalar>&
operator-=(std::vector<Scalar>& x, NoAliasProxy<U> const& y)
{
   subtract(x, y.value());
   return x;
}

} // namespace LinearAlgebra

#endif
