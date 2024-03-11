// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/variant.h
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
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
/* -*- C++ -*- $Id$

  Seralization of boost::variant
*/

#if !defined(VARIANT_H_JSDCHWUI4357784Y7WEHOLWEHO)
#define VARIANT_H_JSDCHWUI4357784Y7WEHOLWEHO

#include "pstream.h"
#include <boost/variant.hpp>
#include <boost/mpl/front.hpp>
#include <boost/mpl/pop_front.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/empty.hpp>

namespace PStream
{

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
opstreambuf<format>& operator<<(opstreambuf<format>& out, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) > const& x);

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstream& operator>>(ipstream& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

// Helper function to load a variant with a given type, specified by the Which parameter.
template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstreambuf<format>& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstream& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

//
// implementation
//

template <int format>
struct variant_save : public boost::static_visitor<>
{
   variant_save(opstreambuf<format>& out_) : out(out_) {}

   template <typename T>
   void operator()(T const& Value) const
   {
      out << Value;
   }
   private:
      opstreambuf<format>& out;
};

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
opstreambuf<format>& operator<<(opstreambuf<format>& out, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) > const& x)
{
   int Which = x.which();
   out << Which;
   boost::apply_visitor(variant_save<format>(out), x);
   return out;
}

template <typename S>
struct variant_load
{
   struct load_null
   {
      template <int format, typename T>
      static void invoke(ipstreambuf<format>&, int, T&)
      {
      }
      template <typename T>
      static void invoke(ipstream&, int, T&)
      {
      }
   };

   struct load_next
   {
      template <int format, typename T>
      static void invoke(ipstreambuf<format>& in, int Which, T& x)
      {
         if (Which == 0)
         {
            typedef typename boost::mpl::front<S>::type head_type;
            head_type Value;
            in >> Value;
            x = Value;
         }
         else
            variant_load<typename boost::mpl::pop_front<S>::type>::load(in, Which-1, x);
      }
      template <typename T>
      static void invoke(ipstream&in, int Which, T& x)
      {
         if (Which == 0)
         {
            typedef typename boost::mpl::front<S>::type head_type;
            head_type Value;
            in >> Value;
            x = Value;
         }
         else
            variant_load<typename boost::mpl::pop_front<S>::type>::load(in, Which-1, x);
      }
   };

   template <int format, typename T>
   static void load(ipstreambuf<format>& in, int Which, T& x)
   {
      typedef typename boost::mpl::eval_if<boost::mpl::empty<S>,
         boost::mpl::identity<load_null>, boost::mpl::identity<load_next> >::type InvokerType;
      InvokerType::invoke(in, Which, x);
   }

   template <typename T>
   static void load(ipstream& in, int Which, T& x)
   {
      typedef typename boost::mpl::eval_if<boost::mpl::empty<S>,
         boost::mpl::identity<load_null>, boost::mpl::identity<load_next> >::type InvokerType;
      InvokerType::invoke(in, Which, x);
   }
};

// In C++17, std::variant has std::variant_size and std::variant_size_v, which seem to have no boost equivalent
template <typename T>
struct variant_size {};

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
struct variant_size<boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>> : public std::integral_constant<std::size_t, boost::mpl::size<typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>::types>::value>
{};

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x)
{
   typedef typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>::types types;
   int Which;
   in >> Which;
   if (Which >= boost::mpl::size<types>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<types>());
   }
   variant_load<types>::load(in, Which, x);
   return in;
}

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstream& operator>>(ipstream& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x)
{
   typedef typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>::types types;
   int Which;
   in >> Which;
   if (Which >= boost::mpl::size<types>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<types>());
   }
   variant_load<types>::load(in, Which, x);
   return in;
}

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstreambuf<format>& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x)
{
   typedef typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>::types types;
   if (Which >= boost::mpl::size<types>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<types>());
   }
   variant_load<types>::load(in, Which, x);
}

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstream& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x)
{
   typedef typename boost::variant<BOOST_VARIANT_ENUM_PARAMS(T)>::types types;
   if (Which >= boost::mpl::size<types>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<types>());
   }
   variant_load<types>::load(in, Which, x);
}

} // namespace PStream

#endif
