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
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
/* -*- C++ -*- $Id$

  Serialization of variant types
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
#include <type_traits>
#include <utility>
#include <variant>

namespace PStream
{

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
opstreambuf<format>& operator<<(opstreambuf<format>& out, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) > const& x);

template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
ipstream& operator>>(ipstream& in, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <int format, typename... T>
opstreambuf<format>& operator<<(opstreambuf<format>& out, std::variant<T...> const& x);

template <int format, typename... T>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, std::variant<T...>& x);

template <typename... T>
ipstream& operator>>(ipstream& in, std::variant<T...>& x);

// Helper function to load a variant with a given type, specified by the Which parameter.
template <int format, BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstreambuf<format>& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
void load_variant(ipstream& in, int Which, boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) >& x);

template <int format, typename... T>
void load_variant(ipstreambuf<format>& in, int Which, std::variant<T...>& x);

template <typename... T>
void load_variant(ipstream& in, int Which, std::variant<T...>& x);

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

template <typename... T>
struct variant_size<std::variant<T...>> : public std::integral_constant<std::size_t, sizeof...(T)>
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

template <std::size_t I, typename Variant, bool Done = (I == std::variant_size<Variant>::value)>
struct std_variant_load;

template <std::size_t I, typename Variant>
struct std_variant_load<I, Variant, false>
{
   template <int format>
   static void load(ipstreambuf<format>& in, int Which, Variant& x)
   {
      if (Which == static_cast<int>(I))
      {
         using value_type = std::variant_alternative_t<I, Variant>;
         value_type Value;
         in >> Value;
         x.template emplace<I>(std::move(Value));
      }
      else
         std_variant_load<I+1, Variant>::load(in, Which, x);
   }

   static void load(ipstream& in, int Which, Variant& x)
   {
      if (Which == static_cast<int>(I))
      {
         using value_type = std::variant_alternative_t<I, Variant>;
         value_type Value;
         in >> Value;
         x.template emplace<I>(std::move(Value));
      }
      else
         std_variant_load<I+1, Variant>::load(in, Which, x);
   }
};

template <std::size_t I, typename Variant>
struct std_variant_load<I, Variant, true>
{
   template <int format>
   static void load(ipstreambuf<format>&, int, Variant&)
   {
   }

   static void load(ipstream&, int, Variant&)
   {
   }
};

template <int format, typename... T>
opstreambuf<format>& operator<<(opstreambuf<format>& out, std::variant<T...> const& x)
{
   if (x.valueless_by_exception())
   {
      PANIC("Cannot serialize valueless std::variant")(tracer::typeid_name<std::variant<T...>>());
   }
   out << static_cast<int>(x.index());
   std::visit([&out](auto const& Value) { out << Value; }, x);
   return out;
}

template <int format, typename... T>
ipstreambuf<format>& operator>>(ipstreambuf<format>& in, std::variant<T...>& x)
{
   int Which;
   in >> Which;
   load_variant(in, Which, x);
   return in;
}

template <typename... T>
ipstream& operator>>(ipstream& in, std::variant<T...>& x)
{
   int Which;
   in >> Which;
   load_variant(in, Which, x);
   return in;
}

template <int format, typename... T>
void load_variant(ipstreambuf<format>& in, int Which, std::variant<T...>& x)
{
   using variant_type = std::variant<T...>;
   if (Which < 0 || static_cast<std::size_t>(Which) >= std::variant_size<variant_type>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<variant_type>());
   }
   std_variant_load<0, variant_type>::load(in, Which, x);
}

template <typename... T>
void load_variant(ipstream& in, int Which, std::variant<T...>& x)
{
   using variant_type = std::variant<T...>;
   if (Which < 0 || static_cast<std::size_t>(Which) >= std::variant_size<variant_type>::value)
   {
      PANIC("Variant out of bounds for type list")(Which)(tracer::typeid_name<variant_type>());
   }
   std_variant_load<0, variant_type>::load(in, Which, x);
}

} // namespace PStream

#endif
