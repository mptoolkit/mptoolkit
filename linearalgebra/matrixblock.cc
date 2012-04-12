// -*- C++ -*- $Id$

namespace LinearAlgebra
{

template <class Scalar, typename Orientation, typername Header, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixBlockRef<Scalar, Orientation, Header, Derived>&>::type
MatrixBlockRef<Scalar, Orientation, Header, Derived>::operator=(U const& x)
{
   // TODO: need a better temporary type
   Matrix<Scalar, Orientation> Temp(x);
   assign(*this, Temp);
   return *this;
}

template <class Scalar, typename Orientation, typename Header, typename Derived>
template <typename U>
inline
typename boost::enable_if<is_matrix<U>, MatrixBlockRef<Scalar, Orientation, Header, Derived>&>::type
MatrixBlockRef<Scalar, Orientation, Header, Derived>::operator=(NoAliasProxy<U> const& x)
{
   assign(*this, x.value());
   return *this;
}

} // namespace LinearAlgebra
