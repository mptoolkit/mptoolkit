/* -*- C++ -*- $Id$

  mapvector.h

  An ordered sparse vector based on std::map

  Created 2005-01-06 Ian McCulloch

  We could improve on the implementation; std::map requires that
  the scalar type is default-constructible, but we don't need that
  for most operations (although others would fail at runtime).
  For example, the set_element, add_element, subtract_element etc interface
  could be used with a scalar type that is not default constructible.
*/

#if !defined(MAPVECTOR_H_JKCHWUITY438T97UY89PU908UWP)
#define MAPVECTOR_H_JKCHWUITY438T97UY89PU908UWP

#include "vectoriteratorfrompair.h"
#include "vectoroperations.h"
#include "vectoraddition.h"
#include "vectortransform.h"
#include <map>
#include "crtp_vector.h"

namespace LinearAlgebra
{

template <typename T>
class MapVector : public VectorBase<MapVector<T> >
{
   private:
      typedef std::map<size_type, T> container_type;
      typedef typename container_type::iterator base_iterator;
      typedef typename container_type::const_iterator base_const_iterator;

   public:
      typedef T value_type;
      typedef VectorIteratorFromPair<base_iterator, vector_iterator_ordered> iterator;
      typedef VectorConstIteratorFromPair<base_const_iterator, 
					  vector_iterator_ordered> const_iterator;
      typedef T& reference;
      typedef T const& const_reference;

      MapVector() : Size_(0) {}

      explicit MapVector(size_type Size) : Size_(Size) {}

      template <typename U>
      MapVector(U const& x, typename boost::enable_if<is_vector<U> >::type* dummy = 0);

      template <typename U>
      MapVector(NoAliasProxy<U> const& x, 
		typename boost::enable_if<is_vector<U> >::type* dummy = 0);

      template <typename U>
      typename boost::enable_if<is_vector<U>, MapVector<T>&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_vector<U>, MapVector<T>&>::type
      operator=(NoAliasProxy<U> const& x);

      const_reference operator[](size_type n) const;
      reference operator[](size_type n);

      size_type size() const { return Size_; }

      const_iterator iterate() const { return const_iterator(Data_.begin(), Data_.end()); }
      const_iterator c_iterate() const { return const_iterator(Data_.begin(), Data_.end()); }
      iterator iterate() { return iterator(Data_.begin(), Data_.end()); }

      iterator iterate_at(size_type n) 
      { return iterator(Data_.find(n), Data_.end()); }

      const_iterator iterate_at(size_type n) const 
      { return const_iterator(Data_.find(n), Data_.end()); }

      size_type nnz() const { return Data_.size(); }

      bool is_zero() const { return Data_.empty(); }

      template <typename U>
      void set_element(size_type n, U const& x);

      template <typename U>
      void add_element(size_type n, U const& x);

      template <typename U>
      void subtract_element(size_type n, U const& x);

      void zero_element(size_type n);

      void zero_all() { Data_.clear(); }

      void resize(size_type NewSize) { Data_.clear(); Size_ = NewSize; }

      void swap(MapVector& Other)
      { std::swap(Size_, Other.Size_); Data_.swap(Other.Data_); }

   private:
      size_type Size_;
      container_type Data_;
};

// interface

template <typename T>
struct interface<MapVector<T> >
{
   //   typedef SEARCH_VECTOR(T, MapVector<T>) type;
   typedef SEARCH_VECTOR(T, void) type;
   typedef T value_type;
};

// iterators

template <typename T>
struct Iterate<MapVector<T>&>
{
   typedef MapVector<T>& argument_type;
   typedef typename MapVector<T>::iterator result_type;
   result_type operator()(argument_type x) const
   {
      return x.iterate();
   }
};

template <typename T>
struct Iterate<MapVector<T> >
{
   typedef MapVector<T> const& argument_type;
   typedef typename MapVector<T>::const_iterator result_type;
   result_type operator()(argument_type x) const
   {
      return x.iterate();
   }
};

// iterate_at

template <typename T>
struct VectorIterateAt<MapVector<T>&>
{
   typedef typename MapVector<T>::iterator result_type;
   typedef MapVector<T>& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type i) const
   {
      return x.iterate_at(i);
   }
};

template <typename T>
struct VectorIterateAt<MapVector<T> >
{
   typedef typename MapVector<T>::const_iterator result_type;
   typedef MapVector<T> const& first_argument_type;
   typedef size_type second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type i) const
   {
      return x.iterate_at(i);
   }
};

// resize

template <typename T>
struct Resize<MapVector<T>&>
{
   typedef void result_type;
   typedef MapVector<T>& first_argument_type;
   typedef size_type second_argument_type;
   void operator()(MapVector<T>& v, size_type n) const
   {
      v.resize(n);
   }
};

// nnz

template <typename T>
struct NNZ<MapVector<T> >
{
   typedef MapVector<T> const& argument_type;
   typedef size_type result_type;

   result_type operator()(argument_type x) const
   {
      return x.nnz();
   }
};

// zero_element - not needed, the INJECTIVE_VECTOR version delegates to the member function

// is_zero

template <typename T>
struct IsZero<MapVector<T> >
{
   typedef MapVector<T> const& argument_type;
   typedef size_type result_type;

   result_type operator()(argument_type x) const
   {
      return x.is_zero();
   }
};

// operator+=, operator-=

template <typename Scalar, typename U>
inline
MapVector<Scalar>&
operator+=(MapVector<Scalar>& x, U const& y)
{
   // TODO: the choice of MapVector as temporary type isn't necessarily optimal here
   MapVector<Scalar> Temp(y);
   add(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
MapVector<Scalar>&
operator+=(MapVector<Scalar>& x, NoAliasProxy<U> const& y)
{
   add(x, y.value());
   return x;
}

template <typename Scalar, typename U>
inline
MapVector<Scalar>&
operator-=(MapVector<Scalar>& x, U const& y)
{
   // TODO: the choice of MapVector as temporary type isn't necessarily optimal here
   MapVector<Scalar> Temp(y);
   subtract(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
MapVector<Scalar>&
operator-=(MapVector<Scalar>& x, NoAliasProxy<U> const& y)
{
   subtract(x, y.value());
   return x;
}

} // namespace LinearAlgebra

#include "mapvector.cc"

#endif
