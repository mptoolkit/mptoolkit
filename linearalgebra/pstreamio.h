/* -*- C++ -*- $Id$

  pstreamio.h

  I/O for vectors and matrices using the PStream framework.

  Created 2005-04-06 Ian McCulloch
*/ 

#if !defined(PSTREAMIO_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define PSTREAMIO_H_DHJFKHURWEIRY4572893475489YRUI34897

namespace LinearAlgebra
{

template <typename Stream, typename T, typename Ti = typename interface<T>::type>
struct SerializeOutInterface {};

template <typename Stream, typename T>
struct SerializeOut : SerializeOutInterface<Stream, T> {};

template <int Format, typename Mat>
typename SerializeOut<PStream::opstreambuf<Format>&, Mat>::result_type
operator<<(PStream::opstreambuf<Format>& out, Mat const& M)
{
   return SerializeOut<PStream::opstreambuf<Format>&, Mat>()(out, M);
}

template <typename Stream, typename T, typename Ti = typename interface<T>::type>
struct SerializeInInterface {};

template <typename Stream, typename T>
struct SerializeIn : SerializeInInterface<Stream, T> {};

template <int Format, typename Mat>
typename SerializeIn<PStream::ipstreambuf<Format>&, Mat&>::result_type
operator>>(PStream::ipstreambuf<Format>& in, Mat& M)
{
   return SerializeIn<PStream::ipstreambuf<Format>&, Mat&>()(in, M);
}

template <int Format, typename Mat>
typename boost::enable_if<is_mutable_proxy<Mat>,
                          SerializeIn<PStream::ipstreambuf<Format>&, Mat&> >::type::result_type
operator>>(PStream::ipstreambuf<Format>& in, Mat const& M)
{
   return SerializeIn<PStream::ipstreambuf<Format>&, Mat&>()(in, const_cast<Mat&>(M));
}

} // namespace LinearAlgebra

#include "pstreamio_vector.h"
#include "pstreamio_matrix.h"

#endif
