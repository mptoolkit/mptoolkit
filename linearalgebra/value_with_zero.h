/* -*- C++ -*- $Id$

  A simple value_with_zero class, for adding a 'universal zero' for types
  that do not support it.

  Created 2005-02-23 Ian McCulloch

  It is not intended that this is a first-class object with respect to
  linear algebra operations, but instead purely as a return type for operators
  that need a result_type that is convertible to T, but also need a zero value.
*/

#if !defined(VALUE_WITH_ZERO_H_SCJKHUTY3879Y3879YT894037OY)
#define VALUE_WITH_ZERO_H_SCJKHUTY3879Y3879YT894037OY

#include <boost/optional.hpp>
#include <boost/mpl/assert.hpp>
#include "interface.h"

namespace LinearAlgebra
{

template <typename T>
class value_with_zero // : private boost::optional<T>
{
   public:
      BOOST_MPL_ASSERT_NOT((has_zero<T>));

      value_with_zero() {}

      value_with_zero(T const& x) : value(boost::optional<T>(x)) {}

      value_with_zero(value_with_zero const& rhs) 
         : value(rhs.value) {}

      template <typename U>
      value_with_zero(value_with_zero<U> const& rhs) 
         : value(rhs.value) {}

      template <typename Expr>
      value_with_zero(Expr const& expr) : value(expr) {}

      value_with_zero& operator=(T const& v) 
        { value = v; return *this; }

      value_with_zero& operator=(value_with_zero const& rhs)
         { value = rhs.value; return *this; }

      template <typename U> 
      value_with_zero& operator=(value_with_zero<U> const& rhs)
          { value = rhs.value; return *this; }

      template <typename Expr>
	 value_with_zero& operator=(Expr const& expr)
	 { value = expr; return *this; }

      T& get() { return value.get(); }
      T const& get() const { return value.get(); }

      bool is_initialized() const { return value.is_initialized(); }

      operator T const&() const { DEBUG_PRECONDITION(!this->is_zero()); 
      DEBUG_PRECONDITION(this->is_initialized()); return this->get(); }
      operator T&() { return this->get(); }

   //      bool is_zero() const { return !static_cast<boost::optional<T> const&>(*this); }
      bool is_zero() const { return !this->is_initialized(); }

      value_with_zero& operator+=(T const& x)
      {
	 if (this->is_zero()) this->operator=(x);
	 else this->get() += x;
	 return *this;
      }

      value_with_zero& operator+=(value_with_zero const& x)
      {
	 if (!x.is_zero())
	 {
	    (*this) += x.get();
	 }
	 return *this;
      }

      template <typename U>
      value_with_zero& operator+=(U const& x)
      {
	 if (this->is_zero()) this->operator=(x);
	 else this->get() += x;
	 return *this;
      }

      template <typename U>
      value_with_zero& operator+=(value_with_zero<U> const& x)
      {
	 if (!x.is_zero())
	 {
	    (*this) += x.get();
	 }
	 return *this;
      }

      value_with_zero& operator-=(T const& x)
      {
	 if (this->is_zero()) this->operator=(-x);
	 else this->get() -= x;
	 return *this;
      }

      value_with_zero& operator-=(value_with_zero const& x)
      {
	 if (!x.is_zero())
	 {
	    (*this) -= x.get();
	 }
	 return *this;
      }

      template <typename U>
      value_with_zero& operator-=(U const& x)
      {
	 if (this->is_zero()) this->operator=(-x);
	 else this->get() -= x;
	 return *this;
      }

      template <typename U>
      value_with_zero& operator-=(value_with_zero<U> const& x)
      {
	 if (!x.is_zero())
	 {
	    (*this) -= x.get();
	 }
	 return *this;
      }

   private:
      boost::optional<T> value;

   template <typename U> friend class value_with_zero;
};

template <typename T>
struct Zero<value_with_zero<T> >
{
   typedef value_with_zero<T> result_type;
   result_type operator()() const { return value_with_zero<T>(); }
};

template <typename T>
struct IsZero<value_with_zero<T> >
{
   typedef bool result_type;
   typedef value_with_zero<T> const& argument_type;
   bool operator()(argument_type x) const
   {
      return x.is_zero();
   }
};

template <typename T, typename Enable>
struct make_value_with_zero
{
   typedef value_with_zero<typename make_value<T>::type> type;
};

template <typename T>
struct make_value_with_zero<
   T
 , typename boost::enable_if<has_zero<typename make_value<T>::type> >::type>
: make_value<T> {};

template <typename T>
struct make_value<value_with_zero<T> > : make_value_with_zero<T> {};

} // namespace LinearAlgebra

#endif
