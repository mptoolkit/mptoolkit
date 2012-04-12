// -*- C++ -*- $Id$

#include "pstreamfwd.h"

namespace Private
{

using PStream::byte_type;

template <EndianType SourceEndian, EndianType DestEndian>
struct EndianHelper;

// If the endianness matches, then pack/unpack is trivial.
template <EndianType E>
struct EndianHelper<E, E>
{
   template <typename Buffer, typename T>
   static void pack(Buffer& B, T const& x)
   {
      B.put_n(reinterpret_cast<byte_type const*>(&x), sizeof(x));
   }

   template <typename Buffer, typename T>
   static void unpack(Buffer& B, T& x)
   {
      B.get_n(reinterpret_cast<byte_type*>(&x), sizeof(T));
   }
};

// If the endian is opposite, then pack/unpack is a simple byte-reverse.
template <>
struct EndianHelper<BigEndian, LittleEndian>
{
   template <typename Buffer, typename T>
   inline static void pack(Buffer& B, T const& x)
   {
      B.template put_bswap<sizeof(T)>(reinterpret_cast<byte_type const*>(&x));
   }

   template <typename Buffer, typename T>
   inline static void unpack(Buffer& B, T& x)
   {
      B.get_n_reverse(reinterpret_cast<byte_type*>(&x), sizeof(T));
   }
};

template <>
struct EndianHelper<LittleEndian, BigEndian> : public EndianHelper<BigEndian, LittleEndian>
{
};

} // namespace private

template <EndianType Endian, typename Buffer, typename T>
inline
void pack_endian(Buffer& B, T const& x)
{
   Private::EndianHelper<CurrentFormat::Endianness, Endian>::pack(B, x);
}

template <EndianType Endian, typename Buffer, typename T>
inline
void unpack_endian(Buffer& B, T& x)
{
   Private::EndianHelper<CurrentFormat::Endianness, Endian>::unpack(B, x);
}
