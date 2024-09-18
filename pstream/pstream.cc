// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/pstream.cc
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "pstreamformats.h"
#include <iterator>
#include <cstring>

namespace PStream
{

namespace Private
{

// The PackList<T> holds, for each type T, an array of pack & unpack functions for each concrete stream type.
template <typename T>
struct PackList
{
   typedef void (*PackFunc)(opstream&, T const&);
   typedef void (*UnpackFunc)(ipstream&, T&);

   static PackFunc PackFuncList[NumConcreteStreamClasses];
   static UnpackFunc UnpackFuncList[NumConcreteStreamClasses];
};

//
// If no streaming function is defined at all, then ordinarily we would
// do a conversion from opstream to opstreambuf<Format> which would then
// implicitly be converted back again, resulting in a recursive call
// and a runtime stack overflow.  Much preferable is a compile-time
// error.  This is achieved by checking the type of the expression
// static_cast<opstreambuf<Format>&>(*Buf.get_buffer()) << Obj
// which, if operator<<(opstreambuf<Format>&, T) is defined, will
// have type opstreambuf<Format>, but if no such function exists
// this will call operator<<(opstream&, T) for a return type of opstream.
// By passing the result to a no-operation template function which
// is defined only for opstreambuf<Format>, a missing inserter will
// result in a compile-time error.
//
// If you see an error "no matching call to
//    templated_opstreambuf_inserter_does_not_exist_for_this_type(opstream)"
// then the stream inserter for the given type is missing.
//

template <int Format>
inline
void templated_opstreambuf_inserter_does_not_exist_for_this_type(opstreambuf<Format>& Stream)
{
}

template <int Format, typename T>
inline
void DoPack(opstream& Buf, T const& Obj)
{
   templated_opstreambuf_inserter_does_not_exist_for_this_type(
      static_cast<opstreambuf<Format>&>(*Buf.get_buffer()) << Obj);
}

//
// the same trick for the extractor that we used for the inserter
//

template <int Format>
inline
void templated_ipstreambuf_extractor_does_not_exist_for_this_type(ipstreambuf<Format>& Stream)
{
}

template <int Format, typename T>
inline
void DoUnpack(ipstream& Buf, T& Obj)
{
   templated_ipstreambuf_extractor_does_not_exist_for_this_type(
   static_cast<ipstreambuf<Format>&>(*Buf.get_buffer()) >> Obj);
}

template <typename T>
typename PackList<T>::PackFunc PackList<T>::PackFuncList[NumConcreteStreamClasses] =
  { DoPack<0>, DoPack<1>, DoPack<2>, DoPack<3> };

template <typename T>
typename PackList<T>::UnpackFunc PackList<T>::UnpackFuncList[NumConcreteStreamClasses] =
  { DoUnpack<0>, DoUnpack<1>, DoUnpack<2>, DoUnpack<3> };

} // namespace Private



namespace Private
{

template <int Format, typename Iter, bool IsTrivialConversion, bool IsSameEndian>
struct WriteHelper  // default case handles non-trivial conversions
{
   typedef typename std::iterator_traits<Iter>::value_type value_type;
   static void write_n(opstreambuf<Format>& Buf, Iter x, size_t N)
   {
      copy_n(x, N, opstreambuf_iterator<Format, value_type>(Buf));
   }
};

// this case handles the case where the conversion is trivial
// and the endianness also matches.  This is just a simple put into the buffer.
// If the iterator is a pointer type, we can further optimize down to a single put_n.
template <int Format, typename Iter>
struct WriteHelper<Format, Iter, true, true>
{
   typedef typename std::iterator_traits<Iter>::value_type value_type;
   static void write_n(opstreambuf<Format>& Buf, Iter x, size_t N)
   {
      for ( ; N > 0; --N)
      {
         Buf.put_n(reinterpret_cast<byte_type const*>(&*x), sizeof(value_type));
         ++x;
      }
   }
};

// optimized case where the iterator type is a pointer
template <int Format, typename T>
struct WriteHelper<Format, T*, true, true>
{
   static void write_n(opstreambuf<Format>& Buf, T* x, size_t N)
   {
      Buf.put_n(reinterpret_cast<byte_type const*>(x), N * sizeof(T));
   }
};

// this case handles the case where the conversion is trivial
// but the endianness is opposite.
template <int Format, typename Iter>
struct WriteHelper<Format, Iter, true, false>
{
   typedef typename std::iterator_traits<Iter>::value_type value_type;
   static void write_n(opstreambuf<Format>& Buf, Iter x, size_t N)
   {
      for (; N > 0; N--)
      {
         Buf.template put_n_bswap<sizeof(value_type)>(reinterpret_cast<byte_type const*>(x), 1);
         ++x;
      }
   }
};

#if 0
template <int Format, typename T>
struct WriteHelper<Format, T const*, true, false>
{
   static void write_n(opstreambuf<Format>& Buf, T const* x, size_t N)
   {
      Buf.template put_n_bswap<sizeof(T)>(reinterpret_cast<byte_type const*>(x), N);
   }
};
#endif

template <int Format, typename T>
struct WriteHelper<Format, T*, true, false>
{
   static void write_n(opstreambuf<Format>& Buf, T* x, size_t N)
   {
      Buf.template put_n_bswap<sizeof(T)>(reinterpret_cast<byte_type const*>(x), N);
   }
};

} // namespace private





//
// implementation of the top-level streaming functions
//

template <typename T>
inline
opstream& operator<<(opstream& Buf, T const& Obj)
{
   Private::PackList<T>::PackFuncList[Buf.format()](Buf, Obj);
   return Buf;
}

template <typename T, int N>
inline
opstream& operator<<(opstream& Buf, T const Obj[N])
{
   Private::PackList<T>::PackFuncList[Buf.format()](Buf, static_cast<T const*>(Obj));
   return Buf;
}

template <typename T>
inline
ipstream& operator>>(ipstream& Buf, T& Obj)
{
   Private::PackList<T>::UnpackFuncList[Buf.format()](Buf, Obj);
   return Buf;
}

//
// opstream
//

inline
opstream::opstream(int Format_, generic_opstreambuf* Buffer_)
  : BaseFormat(Format_), Buffer(Buffer_)
{
}

inline
opstream::opstream(int Format_, byte_type* BufBegin, byte_type* BufEnd, byte_type* BufPtr)
  : BaseFormat(Format_),
    Buffer(generic_opstreambuf::MakeOpstreambuf(Format_, this, BufBegin, BufEnd, BufPtr))
{
}

template <typename T>
inline
void opstream::write(T const& x)
{
   (*this) << x;
}

inline
byte_type* opstream::buf_begin() const
{
   return Buffer->buf_begin();
}

inline
byte_type* opstream::buf_end() const
{
   return Buffer->buf_end();
}

inline
byte_type* opstream::buf_ptr() const
{
   return Buffer->buf_ptr();
}

inline
void opstream::set_buf_ptr(byte_type* Ptr)
{
   Buffer->set_buf_ptr(Ptr);
}

inline
void opstream::set_buffer(byte_type* BufBegin, byte_type* BufEnd, byte_type* BufPtr)
{
   Buffer->set_buffer(BufBegin, BufEnd, BufPtr);
}

inline
void opstream::put_format()
{
   Buffer->put(byte_type(BaseFormat));
}

inline
void opstream::put_format(int NewFormat)
{
   this->set_format(NewFormat);
   Buffer->put(byte_type(NewFormat));
}

//
// ipstream
//

inline
ipstream::ipstream(int Format_, generic_ipstreambuf* Buffer_)
  : BaseFormat(Format_), Buffer(Buffer_)
{
}

inline
ipstream::ipstream(int Format_, byte_type const* BufBegin,
                   byte_type const* BufEnd, byte_type const* BufPtr)
  : BaseFormat(Format_), Buffer(generic_ipstreambuf::MakeIpstreambuf(Format_, this, BufBegin, BufEnd, BufPtr))
{
}

template <typename T>
inline
T ipstream::read()
{
   T temp;
   (*this) >> temp;
   return temp;
}

inline
byte_type const* ipstream::buf_begin() const
{
   return Buffer->buf_begin();
}

inline
byte_type const* ipstream::buf_end() const
{
   return Buffer->buf_end();
}

inline
byte_type const* ipstream::buf_ptr() const
{
   return Buffer->buf_ptr();
}

inline
void ipstream::set_buf_ptr(byte_type const* Ptr)
{
   Buffer->set_buf_ptr(Ptr);
}

inline
void ipstream::set_buf_end(byte_type const* End)
{
   Buffer->set_buf_end(End);
}

inline
void ipstream::set_buffer(byte_type const* BufBegin,
                          byte_type const* BufEnd,
                          byte_type const* BufPtr)
{
   Buffer->set_buffer(BufBegin, BufEnd, BufPtr);
}

inline
int ipstream::get_format()
{
   byte_type NewFormat = Buffer->get();
   this->set_format(NewFormat);
   return NewFormat;
}

//
// generic_opstreambuf
//

inline
void
generic_opstreambuf::set_buffer(byte_type* Begin, byte_type* End, byte_type* Ptr)
{
   BufBegin = Begin;
   BufEnd = End;
   BufPtr = Ptr;
}

inline
void generic_opstreambuf::overflow()
{
   BackPointer->overflow();
}

inline
void generic_opstreambuf::flush()
{
   BackPointer->flush();
}

inline
void generic_opstreambuf::put(byte_type b)
{
   while (BufPtr == BufEnd) this->overflow();
   *BufPtr++ = b;
}

inline
void generic_opstreambuf::put_n(byte_type const* Src, size_t Count)
{
   size_t Avail = this->buf_avail();
   while (Avail < Count)
   {
      std::memcpy(BufPtr, Src, Avail);
      BufPtr = BufEnd;             // same as BufPtr += Avail
      Src += Avail;
      Count -= Avail;
      this->overflow();
      Avail = this->buf_avail();
   }
   std::memcpy(BufPtr, Src, Count);
   BufPtr += Count;
}

/*
inline
void generic_opstreambuf::put_n_reverse(byte_type const* Src, size_t Count)
{
   // simple-minded implementation
   byte_type const* Ptr = Src + Count;
   while (Ptr != Src)
   {
      --Ptr;
      this->put(*Ptr);
   }
}
*/

template <int Atom>
inline
void generic_opstreambuf::put_bswap(byte_type const* Src)
{
   Src += Atom;
   for (int i = 0; i < Atom; ++i)
   {
      this->put(*--Src);
   }
}

template <int Atom>
inline
void generic_opstreambuf::put_bswap_unchecked(byte_type const* Src)
{
   DEBUG_PRECONDITION(this->buf_avail() >= Atom);
   Src += Atom;
   for (int i = 0; i < Atom; ++i)
   {
      *BufPtr++ = *--Src;
   }
}

template <int Atom>
inline
void generic_opstreambuf::put_n_bswap(byte_type const* Src, size_t Count)
{
   while (Count > 0)
   {
      size_t Avail = this->buf_avail();

      if (Avail == 0)
      {
         this->overflow();
         Avail = this->buf_avail();
      }

      if (Avail > Count * Atom) Avail = Count * Atom;
      while (Avail >= Atom)
      {
         this->put_bswap_unchecked<Atom>(Src);
         Src += Atom;
         Avail -= Atom;
         --Count;
      }

      DEBUG_CHECK(Avail < Atom);
      DEBUG_CHECK(Avail == 0 || Count > 0);

      if (Avail > 0)
      {
         this->put_bswap<Atom>(Src);
         Src += Atom;
         --Count;
      }
   }
}

inline
void generic_opstreambuf::bump(size_t Count)
{
   DEBUG_CHECK(Count <= this->buf_avail());
   BufPtr += Count;
}

//
// opstreambuf
//

template <int Format>
opstreambuf<Format>::opstreambuf(opstream* BackPointer_,
                                byte_type* Beg, byte_type* End, byte_type* Ptr)
  : generic_opstreambuf(BackPointer_, Beg, End, Ptr)
{
}

template <int Format>
opstreambuf<Format>::~opstreambuf()
{
}

template <int Format>
template <typename T>
inline
void opstreambuf<Format>::write(T const& x)
{
   (*this) << x;
}

template <int Format, typename T, bool IsFundamental = pstream_type_traits<T>::is_fundamental>
struct ConversionTraitsHelper
{
   static bool const trivial_conversion = false;
};

template <int Format, typename T>
struct ConversionTraitsHelper<Format, T, true>
{
   typedef typename format_traits<Format>::Format format;
   typedef typename format::template TypeTraits<T> Traits;
   static bool const trivial_conversion =
      CurrentFormat::IsTrivialConversion<T, typename Traits::type>::value;
};


template <int Format>
template <typename T>
inline
void opstreambuf<Format>::write_n(T const* x, size_t N)
{
   typedef typename format_traits<Format>::Format format;
   Private::WriteHelper<Format, T const*,
     ConversionTraitsHelper<Format, T>::trivial_conversion,
     CurrentFormat::Endianness == format::Endianness>::write_n(*this, x, N);
}

//
// generic_ipstreambuf
//


inline
byte_type generic_ipstreambuf::get()
{
   while (BufPtr == BufEnd) underflow();
   return *BufPtr++;
}

inline
void generic_ipstreambuf::get_n(byte_type* Dest, size_t Count)
{
   size_t Avail = this->buf_avail();
   while (Avail < Count)
   {
      std::memcpy(Dest, BufPtr, Avail);
      BufPtr = BufEnd;             // same as BufPtr += Avail
      Count -= Avail;
      Dest += Avail;
      this->underflow();
      Avail = this->buf_avail();
   }
   std::memcpy(Dest, BufPtr, Count);
   BufPtr += Count;
}

inline
void generic_ipstreambuf::get_n_reverse(byte_type* Dest, size_t Count)
{
   // simple-minded implementation
   byte_type* Ptr = Dest + Count;
   while (Ptr != Dest)
   {
      --Ptr;
      *Ptr = this->get();
   }
}

inline
void generic_ipstreambuf::set_buffer(byte_type const* Begin, byte_type const* End, byte_type const* Ptr)
{
   BufBegin = Begin;
   BufEnd = End;
   BufPtr = Ptr;
}

inline
void generic_ipstreambuf::bump(size_t Count)
{
   DEBUG_CHECK(Count <= this->buf_avail());
   BufPtr += Count;
}

inline
void generic_ipstreambuf::underflow()
{
   BackPointer->underflow();
}

//
// ipstreambuf
//

template <int Format>
ipstreambuf<Format>::ipstreambuf(ipstream* BackPointer,
                         byte_type const* Beg, byte_type const* End, byte_type const* Ptr)
  : generic_ipstreambuf(BackPointer, Beg, End, Ptr)
{
}

template <int Format>
ipstreambuf<Format>::~ipstreambuf()
{
}

template <int Format>
template <typename T>
inline
T ipstreambuf<Format>::read()
{
   T temp;
   (*this) >> temp;
   return temp;
}

//
// versioning
//

inline
VersionSentry::VersionSentry(ipstream& In, VersionTag const& Tag_, int v_)
   : Tag(Tag_), in(&In), out(NULL), v(v_)
{
   In.push_version(Tag, v);
}

template <int Format>
inline
VersionSentry::VersionSentry(ipstreambuf<Format>& In, VersionTag const& Tag_, int v_)
   : Tag(Tag_), in(&In), out(NULL), v(v_)
{
   In.push_version(Tag, v);
}

inline
VersionSentry::VersionSentry(opstream& Out, VersionTag const& Tag_, int v_)
   : Tag(Tag_), in(NULL), out(&Out), v(v_)
{
   Out.push_version(Tag, v);
}

template <int Format>
inline
VersionSentry::VersionSentry(opstreambuf<Format>& Out, VersionTag const& Tag_, int v_)
   : Tag(Tag_), in(NULL), out(&Out), v(v_)
{
   Out.push_version(Tag, v);
}

inline
VersionSentry::~VersionSentry()
{
   if (in)
      in->pop_version(Tag);
   if (out)
      out->pop_version(Tag);
}

inline
void
VersionSentry::change_version(int NewV)
{
   v = NewV;
   if (in)
   {
      in->pop_version(Tag);
      in->push_version(Tag, NewV);
   }
   if (out)
   {
      out->pop_version(Tag);
      out->push_version(Tag, NewV);
   }
}

//
// Polymorphic streaming support
//

template <typename T>
class ClassStreamableFactory : public StreamableFactory
{
   public:
      ClassStreamableFactory() {}

      StreamableBase* construct(ipstream& stream);
};

template <typename T>
StreamableBase* ClassStreamableFactory<T>::construct(ipstream& stream)
{
   StreamableBase* Ret = new T;
   Ret->ReadStream(stream);
   return Ret;
}



template <typename T>
class ClassStreamableViaCtorFactory : public StreamableFactory
{
   public:
      ClassStreamableViaCtorFactory() {}

      StreamableBase* construct(ipstream& stream);
};

template <typename T>
StreamableBase* ClassStreamableViaCtorFactory<T>::construct(ipstream& stream)
{
   return new T(stream);
}

template <typename T>
void Register(const std::string& s)
{
   RegisterClass(typeid(T), s, new ClassStreamableFactory<T>());
}

//
// pre-defined streaming operators for non-fundamental types.
//

//
// copy_n is useful - but it should be defined somewhere else!
//

template <typename InIter, typename OutIter>
OutIter copy_n(InIter In, size_t Count, OutIter Out)
{
   for (size_t s = 0; s < Count; ++s)
   {
      *Out++ = *In++;
   }
   return Out;
}

template <int Format>
opstreambuf<Format>& operator<<(opstreambuf<Format>& out, char const* String)
{
   typename opstreambuf<Format>::size_type Size = std::strlen(String);
   out << Size;
   out.put_n(reinterpret_cast<byte_type const*>(String), Size);
   return out;
}

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::basic_string<T> const& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   using std::copy;
   copy(vec.begin(), vec.end(), opstreambuf_iterator<Format, T>(stream));
   return stream;
}

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::basic_string<T>& vec)
{
   typename ipstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.resize(Size);
   copy_n(ipstreambuf_iterator<Format, T>(stream), Size, vec.begin());
   //   TRACE(Size)(vec);
   return stream;
}

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::vector<T> const& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;

   if (Size > 0)
     stream.write_n(&vec[0], Size);
   return stream;
}

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::vector<T>& vec)
{
   typename ipstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.resize(Size);
   copy_n(ipstreambuf_iterator<Format, T>(stream), Size, &vec[0]);
   return stream;
}

template <int Format>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::vector<bool> const& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   for (int i = 0; i < Size; ++i)
   {
      stream << vec[i];
   }
   return stream;
}

template <int Format>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::vector<bool>& vec)
{
   typename opstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.clear();
   for (int i = 0; i < Size; ++i)
   {
      bool b;
      stream >> b;
      vec.push_back(b);
   }
   return stream;
}

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::deque<T>& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   copy(vec.begin(), vec.end(), opstreambuf_iterator<Format, T>(stream));
   return stream;
}

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::deque<T>& vec)
{
   typename opstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.clear();
   copy_n(ipstreambuf_iterator<Format, T>(stream), Size, std::back_inserter(vec));
   return stream;
}

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::list<T>& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   using std::copy;
   copy(vec.begin(), vec.end(), opstreambuf_iterator<Format, T>(stream));
   return stream;
}

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::list<T>& vec)
{
   typename opstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.clear();
   copy_n(ipstreambuf_iterator<Format, T>(stream), Size, std::back_inserter(vec));
   return stream;
}

template <int Format, typename T, class Comp>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::set<T, Comp> const& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   using std::copy;
   copy(vec.begin(), vec.end(), opstreambuf_iterator<Format, T>(stream));
   return stream;
}

template <int Format, typename T, class Comp>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::set<T, Comp>& vec)
{
   typename opstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.clear();
   copy_n(ipstreambuf_iterator<Format, T>(stream), Size, std::inserter(vec, vec.end()));
   return stream;
}

template <int Format, typename T1, typename T2, class Comp>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::map<T1, T2, Comp>& vec)
{
   typename opstreambuf<Format>::size_type Size = vec.size();
   stream << Size;
   using std::copy;
   copy(vec.begin(), vec.end(), opstreambuf_iterator<Format, std::pair<T1, T2> >(stream));
   return stream;
}

template <int Format, typename T1, typename T2, class Comp>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::map<T1, T2, Comp>& vec)
{
   typename opstreambuf<Format>::size_type Size;
   stream >> Size;
   vec.clear();
   while (Size > 0)
   {
      std::pair<T1, T2> x;
      stream >> x;
      vec.insert(x);
      --Size;
   }
   //   copy_n(ipstreambuf_iterator<Format, std::pair<T1, T2> >(stream), Size, std::inserter(vec));
   return stream;
}

template <int Format, typename T1, typename T2>
inline
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::pair<T1, T2> const& Obj)
{
   stream << Obj.first;
   stream << Obj.second;
   return stream;
}

template <int Format, typename T1, typename T2>
inline
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::pair<T1, T2>& Obj)
{
   stream >> Obj.first;
   stream >> Obj.second;
   return stream;
}

template <int Format, typename T>
inline
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::complex<T> const& Obj)
{
   stream << Obj.real();
   stream << Obj.imag();
   return stream;
}

template <int Format, typename T>
inline
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::complex<T>& Obj)
{
   //  T r,i;
  stream >> reinterpret_cast<T&>(Obj);
  stream >> *(reinterpret_cast<T*>(&Obj)+1);
  //  Obj = std::complex<T>(r,i);
  return stream;
}


} // namespace PStream
