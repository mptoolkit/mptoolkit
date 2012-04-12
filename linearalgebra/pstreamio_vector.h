/* -*- C++ -*- $Id$

  pstreamio_vector.h

  I/O for vectors using the PStream framework.

  Created 2005-05-12 Ian McCulloch
*/ 

#if !defined(PSTREAMIO_VECTOR_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define PSTREAMIO_VECTOR_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "vectoroperations.h"
#include "pstream/pstream.h"

namespace LinearAlgebra
{

// The available formats for saving vectors.
// This isn't an enumeration because we need to control the representation
// on disk.
struct VectorFormats
{
   typedef char type;
   static type const Dense = 'm';
   static type const Coordinate = 'c';
};

template <int Format, typename Vec, typename VecV, typename VecI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Vec,
                             DENSE_VECTOR(VecV, VecI)>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Vec const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type V) const
   {
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Size = size(V);
#if !defined(PSTREAMIO_OLD_FORMAT)
      out << VectorFormats::Dense;
#endif
      out << Size;

      typename const_iterator<Vec>::type I = iterate(V);
      while (I)
      {
         out << *I;
         ++I;
      }
      return out;
   }
};

template <int Format, typename Vec, typename VecV, typename VecI>
struct SerializeOutInterface<PStream::opstreambuf<Format>&,
                             Vec,
                             COMPRESSED_VECTOR(VecV, VecI)>
{
   typedef PStream::opstreambuf<Format>& result_type;
   typedef PStream::opstreambuf<Format>& first_argument_type;
   typedef Vec const& second_argument_type;

   result_type operator()(first_argument_type out, second_argument_type V) const
   {
#if defined(PSTREAMIO_OLD_FORMAT)
      PANIC("Cannot serialize sparse vectors with PSTREAMIO_OLD_FORMAT");
#endif
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Size = size(V), Nnz = nnz(V);
      out << VectorFormats::Coordinate;
      out << Size << Nnz;

      typename const_iterator<Vec>::type I = iterate(V);
      while (I)
      {
         st r = I.index();
         out << r;
         out << *I;
         ++I;
      }
      return out;
   }
};



namespace Private
{

//
// implementation helpers
//

template <typename Stream, typename Vec>
void do_serialize_in_dense_vector(Stream& in, Vec& V)
{
   typename iterator<Vec>::type I = iterate(V);
   while (I)
   {
      in >> *I;
      ++I;
   }
}

template <int Format, typename Vec>
void do_serialize_in_coordinate_vector(PStream::ipstreambuf<Format>& in, Vec& V)
{
   typedef typename PStream::ipstreambuf<Format>::size_type st;
   typedef typename interface<Vec>::value_type value_t;
   st nnz;
   in >> nnz;
   for (st i = 0; i < nnz; ++i)
   {
      st r;
      in >> r;
      add_element(V, r, in.template read<value_t>());
   }
}

template <int Format, typename Vec>
void do_serialize_in_coordinate_vector(PStream::ipstreambuf<Format>& in, Vec const& V)
{
   do__serialize_in_coordinate_vector(in, const_cast<Vec&>(V));
}

} // namespace Private

template <int Format, typename Vec, typename VecV, typename VecI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Vec&,
                            DENSE_VECTOR(VecV, VecI)>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Vec& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type V) const
   {
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Size;
      VectorFormats::type t;
#if defined(PSTREAMIO_OLD_FORMAT)
      t = VectorFormats::Dense;
#else
      in >> t;
#endif
      in >> Size;
      try_resize(V, Size);
      if (t == VectorFormats::Dense)
      {
         Private::do_serialize_in_dense_vector(in, V);
      }
      //      else if (t == VectorFormats::Coordinate)
      //      {
      //         zero_all(V);
      //         Private::do_serialize_in_coordinate_vector(in, V);
      //      }
      else
      {
         PANIC("Unsupported vector format.")(t)(int(t));
      }
      return in;
   }
};

template <int Format, typename Vec, typename VecV, typename VecI>
struct SerializeInInterface<PStream::ipstreambuf<Format>&,
                            Vec&,
                            COMPRESSED_VECTOR(VecV, VecI)>
{
   typedef PStream::ipstreambuf<Format>& result_type;
   typedef PStream::ipstreambuf<Format>& first_argument_type;
   typedef Vec& second_argument_type;

   result_type operator()(first_argument_type in, second_argument_type V) const
   {
#if defined(PSTREAMIO_OLD_FORMAT)
      PANIC("Cannot serialize sparse vectors with PSTREAMIO_OLD_FORMAT");
#endif
      typedef typename PStream::opstreambuf<Format>::size_type st;
      st Size;
      VectorFormats::type t;
      in >> t;
      in >> Size;
      try_resize(V, Size);
      if (t == VectorFormats::Coordinate)
      {
         zero_all(V);
         Private::do_serialize_in_coordinate_vector(in, V);
      }
      else
      {
         PANIC("Unsupported format for sparse vector.")(t)(int(t));
      }
      return in;
   }
};

} // namespace LinearAlgebra

#endif
