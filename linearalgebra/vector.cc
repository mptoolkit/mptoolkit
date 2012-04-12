// -*- C++ -*- $Id$

#include "vectormemproxy.h"

namespace LinearAlgebra
{

template <class Scalar>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, VectorRef<Scalar>&>::type
VectorRef<Scalar>::operator=(U const& x)
{
   Vector<Scalar> Temp(x);
   assign(*this, Temp);
   return *this;
}

template <class Scalar>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, VectorRef<Scalar>&>::type
VectorRef<Scalar>::operator=(NoAliasProxy<U> const& x)
{
   assign(*this, x.value());
   return *this;
}

template <class Scalar>
inline
void 
VectorRef<Scalar>::resize(size_type NewSize)
{
   // Early return in the case that the data block is not shared
   // and we don't need a resize.
   size_type Size = this->size();
   if (NewSize == Size && !Block.is_shared()) return;

   DataBlock<Scalar> NewData(NewSize);
   Block.data() = NewData;
}

template <class Scalar>
inline
Vector<Scalar>::Vector(size_t Size_)
  : VectorRef<Scalar>(Size_)
{
}

template <class Scalar>
inline
Vector<Scalar>::Vector(size_t Size_, Scalar const& fill_)
  : VectorRef<Scalar>(Size_, fill_)
{
}

template <class Scalar>
template <typename Iter>
inline
Vector<Scalar>::Vector(Iter first, Iter last, 
		       typename boost::enable_if<
		       boost::mpl::not_<boost::is_arithmetic<Iter> > >::type*)
   : VectorRef<Scalar>(std::distance(first, last))
{
   assign(*this, std::vector<Scalar>(first, last));
}

template <class Scalar>
template <typename U>
inline
Vector<Scalar>::Vector(U const& x, 
		       typename boost::enable_if<is_vector<U> >::type*)
  : VectorRef<Scalar>(Size<U>()(x))
{
   assign(*this, x);
}

template <class Scalar>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, Vector<Scalar>&>::type
Vector<Scalar>::operator=(U const& x)
{
   // We must use a temporary here, as the right hand side could
   // be aliased without increasing the reference count
   Vector<Scalar> Temp(x);
   //this->resize(Temp.size());
   assign(*this, Temp);
   return *this;
}

template <class Scalar>
template <typename U>
inline
typename boost::enable_if<is_vector<U>, Vector<Scalar>&>::type
Vector<Scalar>::operator=(NoAliasProxy<U> const& x)
{
   this->resize(Size<typename NoAliasProxy<U>::value_type>()(x.value()));
   assign(*this, x.value());
   return *this;
}

template <class Scalar>
inline
Vector<Scalar>&
Vector<Scalar>::operator=(Vector<Scalar> const& V)
{
   this->block() = V.block().copy();
   return *this;
}

#if 0
template <typename Scalar, typename RHS, typename Any, typename RValue, typename U>
struct Assign<Vector<Scalar>, RHS, Any, LOCAL_VECTOR(RValue, U)>
{
   static void apply(Vector<Scalar>& x, RHS const& y)
   {
      fill(x, zero<Scalar>());
      typename const_iterator<RHS>::type Iter = iterate(y);
      while (y)
      {
	 x[y.index()] = *y;
	 ++y;
      }
   }
};

// this overload is redundant...
template <typename Scalar, typename RHS, typename RValue, typename U>
struct Assign<Vector<Scalar>, RHS, CONTIGUOUS_VECTOR(Scalar, Vector<Scalar>), DENSE_VECTOR(RValue, U)>
{
   static void apply(Vector<Scalar>& x, RHS const& y)
   {
      iter_copy(iterate(y), iterate(x));
   }
};
#endif

} // namespace LinearAlgebra
