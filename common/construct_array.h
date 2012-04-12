// -*- C++ -*- $Id$

/*
  construct_array.h

  Functions to construct and destroy arrays at preallocated memory.

  Created 2004-04-19 Ian McCulloch

  Synopsis:
    The array placement-new construct in C++98 is essentially unuseable.
    construct_array() and destroy_array() are intended to be substitites, for
    constructing an array at a fixed memory address.

*/

#if !defined(CONSTRUCT_ARRAY_H_JHU438R9U89UF89YPJJF89PREUHP)
#define CONSTRUCT_ARRAY_H_JHU438R9U89UF89YPJJF89PREUHP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <memory>

namespace ext
{

//
// construct_array
//
// default-constructs an array of type T, with Size elements.  
// buf is assumed to be uninitialized memory of length at least sizeof(T) * Size.  
// If any of the constructors throws an exception, destructors are called
// for all fully constructed elements.
//

template <typename T>
void construct_array(T* buf, size_t Size);

//
// copy-constructs an array of type T, with Size elements, into uninitialized memory.
//

template <typename T>
void construct_array(T* buf, size_t Size, T const& x);

//
// destroy_array
//
// invokes destructors for the array of Size elements located at buf.
//

template <typename T>
void destroy_array(T const* buf, size_t Size);

//
// elide_construction
//
// a metafunction that, if set to true_, leaves the memory of an array
// uninitialized.
//

template <typename T>
struct elide_construction;

//
// elide_destruction
//
// a metafunction that, if set to true_, causes the destructors for
// an array to be not called.
//

template <typename T>
struct elide_destruction;

//
// hacks for std::complex
//

#if defined(NDEBUG)
template <typename T>
struct elide_construction : boost::has_trivial_constructor<T> {};

template <typename T>
struct elide_destruction : boost::has_trivial_destructor<T> {};

template <typename T>
struct elide_construction<std::complex<T> > : boost::has_trivial_constructor<T> {};

template <typename T>
struct elide_destruction<std::complex<T> > : boost::has_trivial_destructor<T> {};
#else

template <typename T>
struct elide_construction : boost::mpl::false_ {};

template <typename T>
struct elide_destruction : boost::mpl::false_ {};
#endif


// implementation


namespace Private
{

template <typename T, typename Enable = void>
struct DestroyHelper
{
   static void apply(T const* buf, size_t Size);
};

template <typename T>
struct DestroyHelper<T, typename boost::enable_if<elide_destruction<T> >::type>
{
   static void apply(T const* buf, size_t Size) { }
};

template <typename T, typename Enable>
void DestroyHelper<T, Enable>::apply(T const* buf, size_t Size)
{
   // destroy in reverse order to construction
   for ( ; Size != 0; --Size)
   {
     buf[Size-1].~T();
   }
}

template <typename T, typename Enable = void>
struct ConstructHelper
{
   static void apply(T* buf, size_t Size);
   static void apply(T* buf, size_t Size, T const& x);
};

inline
double nans(char const*)
{
   unsigned long long x = 0x7ff7ffffffffffffULL;
   char const* c = static_cast<char const*>(static_cast<void const*>(&x));
   return *static_cast<double const*>(static_cast<void const*>(c));
}

template <typename T>
struct DebugInit
{
   static T apply() { return T(); }
};

template <>
struct DebugInit<double>
{
   static double apply() { return nans(""); }
};

template <>
struct DebugInit<std::complex<double> >
{
   static std::complex<double> apply() { return std::complex<double>(nans(""), nans("")); }
};

template <typename T>
struct ConstructHelper<T, typename boost::enable_if<elide_construction<T> >::type>
{
   static void apply(T* buf, size_t Size) { }
   static void apply(T* buf, size_t Size, T const& x)
   { std::uninitialized_fill(buf, buf+Size, x); }
};

#if defined(NDEBUG)
template <typename T, typename Enable>
void ConstructHelper<T, Enable>::apply(T* buf, size_t Size)
{
   size_t i = 0;
   try
   {
      for ( ; i < Size; ++i)
      {
	 new (buf+i) T;
      }
   }
   catch (...)
   {
      // if buf[i]'s constructor threw, then destuct all of the previous objects
      for ( ; i != 0; --i)
      {
	 buf[i-1].~T();
      }
      throw;
   }
}
#else
template <typename T, typename Enable>
void ConstructHelper<T, Enable>::apply(T* buf, size_t Size)
{
   apply(buf, Size, DebugInit<T>::apply());
}
#endif

template <typename T, typename Enable>
void ConstructHelper<T, Enable>::apply(T* buf, size_t Size, T const& x)
{
   size_t i = 0;
   try
   {
      for ( ; i < Size; ++i)
      {
	 new (buf+i) T(x);
      }
   }
   catch (...)
   {
      // if buf[i]'s constructor threw, then destuct all of the previous objects
      for ( ; i != 0; --i)
      {
	 buf[i-1].~T();
      }
      throw;
   }
}
} // namespace Private

template <typename T>
void construct_array(T* buf, size_t Size)
{
   Private::ConstructHelper<T>::apply(buf, Size);
}

template <typename T>
void construct_array(T* buf, size_t Size, T const& x)
{
   Private::ConstructHelper<T>::apply(buf, Size, x);
}

template <typename T>
void destroy_array(T const* buf, size_t Size)
{
   Private::DestroyHelper<T>::apply(buf, Size);
}

} // namespace ext

#endif
