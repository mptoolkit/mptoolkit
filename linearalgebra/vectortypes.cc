// -*- C++ -*- $Id$

#include "dataops.h"

namespace LinearAlgebra
{

template <class Scalar>
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

#if defined(USE_PSTREAM)
template <int Format, typename T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& out, 
					 Vector<T> const& Vec)
{
   typename PStream::opstreambuf<Format>::size_type Size = Vec.size();
   out << Size;
   PStream::copy_n(Vec.begin(), Vec.size(), PStream::opstreambuf_iterator<Format, T>(out));
   return out;
}

template <int Format, typename T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& in, Vector<T>& Vec)
{
   typename PStream::ipstreambuf<Format>::size_type Size;
   in >> Size;
   Vec.resize(Size, false);
   PStream::copy_n(PStream::ipstreambuf_iterator<Format, T>(in), Size, Vec. begin());
   return in;
}
#endif

template <class Scalar>
Vector<Scalar>::Vector(size_t Size_)
  : VectorRef<Scalar>(Size_)
{
}

template <class Scalar>
Vector<Scalar>::Vector(size_t Size_, Scalar const& fill_)
  : VectorRef<Scalar>(Size_, fill_)
{
}

template <class Scalar>
template <typename Iter>
Vector<Scalar>::Vector(Iter first, Iter last, 
		       typename boost::enable_if<
		       boost::mpl::not_<boost::is_arithmetic<Iter> > >::type*)
  : VectorRef<Scalar>(last-first)
{
   // FIXME: this really should call an overload of assign()
   ops::fast_copy(first, last, this->begin());
}

#if 0
template <class Scalar>
template <typename S, typename D>
Vector<Scalar>::Vector(VectorConstExpression<S, D> const& V)
  : VectorRef<Scalar>(V.size())
{
   assign(*this, V);
}
#endif

#if 0
template <class Scalar>
template <class S2, class D2>
Vector<Scalar>&
Vector<Scalar>::operator=(VectorConstExpression<S2 const, D2> const& V)
{
   this->resize(V.size(), false);
   assign(*this, V);
   return *this;
}
#endif

template <class Scalar>
Vector<Scalar>&
Vector<Scalar>::operator=(Vector<Scalar> const& V)
{
   this->block() = V.block().copy();
   return *this;
}

} // namespace LinearAlgebra
