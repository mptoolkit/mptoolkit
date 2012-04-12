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
   };

   template <int format, typename T>
   static void load(ipstreambuf<format>& in, int Which, T& x)
   {
      typedef typename boost::mpl::eval_if<boost::mpl::empty<S>,
         boost::mpl::identity<load_null>, boost::mpl::identity<load_next> >::type InvokerType;
      InvokerType::invoke(in, Which, x);
   }
};

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

} // namespace PStream

#endif
