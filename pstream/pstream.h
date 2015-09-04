/* -*- C++ -*- $Id$

  The new version of pstream.h

  Created 2002-04-08 Ian McCulloch

  Object serialization library.

  This allows different formats differing in representation (endianness, size, or even
  binary versus eg. XML), in a way that is still
  fairly efficient.  The heirachy has similarities to the standard IOStreams library,
  but also some differences.

  The concrete stream classes are derived from opstream/ipstream.  These can be used in
  a way that is exactly analagous to IOStreams.

  opstream/ipstream contain an instance of generic_[i|o]pstreambuf.  The generic_pstreambuf
  class contains pointers to the actual buffer, and has members
  generic_opstreambuf:
    (buffer manipulation)
     byte_type* buf_begin();
     byte_type* buf_end();
     byte_type* buf_ptr();
     size_t buf_avail();
     void set_buffer(byte_type* Begin, byte_type* End, byte_type* Ptr);

    (put)
     void put(byte_type);
     void put_n(byte_type* Src, size_t Count);
     void put_n_reverse(byte_type* Src, size_t Count);

     void overflow();

  generic_ipstreambuf:
    (buffer manipulation)
     byte_type* buf_begin();
     byte_type* buf_end();
     byte_type* buf_ptr();
     size_t buf_avail();
     void set_buffer(byte_type* Begin, byte_type* End, byte_type* Ptr);

    (get)
     byte_type get();
     void get_n(byte_type* Dest, size_t Count);
     void get_n_reverse(byte_type* Dest, size_t Count);

     void underflow();

  The buffer is controlled by the associated stream class.  Overflow() and underflow()
  delegate to the corresponding pure virtual functions in opstreambuf/ipstreambuf.
  To do this, the buffer classes hold a back-pointer to the stream class.

  The format (representation) is determined by the type of the actual buffer class.
  These are labelled by integers:

  template <int Format>
  class opstreambuf : public generic_opstreambuf { ... };

  (The reasoning for choosing an int here is that some versions of g++ seem
  to have problems with integral template parameters of other types.  int at least
  appears safe.)

  For user types, there are two possibilities for writing the streaming operators.
  The easiest is to use the opstream/ipstream classes as usual.  eg,

  opstream& operator<<(opstream& out, MyType const& x) { ... }

  Alternatively, it is possible to write a template function that streams to the buffer directly.
  This is done by writing a template function

  template <int Format>
  opstreambuf<Format>& operator<<(opstreambuf<Format>& out, MyType const& x) { ... }

  The choice as to which function to use depends on the circumstances.  All non-pointer builtin
  types are handled by template functions, so streaming an aggregate of builtins is
  faster if it is itself handled by a template function.  The penalty is more
  code bloat: the streaming functions will be instantiated for each possible data format.
  Also, using the pstreambuf interface to stream an aggregate of objects that individually use the
  pstream interface gives no speedup.

  Streaming operators for user types should not directly access the buffer (that would bypass
  the formatting), but instead simply stream its component types.  Streaming operators
  are pre-defined for all non-pointer builtin types, as well as the fixed format
  [u]int8, [u]int16, [u]int32 and [u]int64 types, defined in common/inttype.h.
  These types always retain the advertized size, irrespective of the platform or format.

  WARNING: It is not possible to portably stream types which are platform-dependent typedef's,
  such as size_t, ptrdiff_t, ssize_t, off_t, int8_t, int16_t, int32_t, int64_t etc.
  Always convert these types to either a known builtin (int, long etc),
  or use one of the int8, int16, int32, int64 etc types defined in inttype.h,
  or use one of the format types defined in [i|o]pstreambuf<Format>::size_type, difference_type etc.

  The int template parameter of the pstreambuf classes is for efficiency: rather than 
  using virtual functions, the pstream class contains the int parameter of its concrete buffer class,
  and dispatch to the correct format object is done by an array lookup.  In pseudo-code:

  template <typename T>
  opstream& operator<<(opstream& out, T const& x)
  {
     typedef opstreambuf<out.format()> concrete_streambuf_type;
     static_cast<concrete_streambuf_type&>(*out.buffer()) << x;
     return out;
  }

  This is not real C++, because out.format() is a runtime variable, but what is needed
  for the opstreambuf is a compile-time template parameter.  This trickery is handled
  by the PStream::Private::PackList structure defined in pstream.cc.  The end result
  is the same though.  The tradeoff is that the number of concrete buffer classes
  (ie. number of different output representations) must be a compile-time constant
  and all of them must be known at the point of instantiation of the streaming operator.
  
  The [i|o]pstreambuf<> classes have an implicit conversion to [i|o]pstream,
  which handles the case where a [i|o]pstreambuf extractor/inserter contains
  an object that has only an [i|o]pstream operator.

  NOTE: be careful when chaining calls inside [i|o]pstreambuf extractor/inserter functions.
  If a call involves a conversion to [i|o]pstream, then all calls to the right of this
  will act on the pstream, and not the pstreambuf!  This is a serious problem that I'm
  unsure how to fix: workaround is to not use chaining of insertion/extraction operators
  on pstreambufs.  Chaining on pstream is fine.

  The data formats are orthogonal to the concrete stream types (ie. the concrete class
  derived from [i|o]pstream that implements the overflow()/underflow() function).
  In principle, every data format is available to every stream type, the
  choice of data formats is controlled by the concrete stream class.  For example,
  a concrete stream class might have a function which changes the data format at
  any time, or a concrete stream class might have a fixed data format.

  Although the number of data formats (number of pstreambuf classes that will be instantiated)
  is fixed at compile time, the number of stream types derived from [i|o]pstream is not fixed at all.

  TODO: Explanation of how streaming of polymorphic objects works.  Short answer: inherit from
  StreamableBase and implement the ReadStream() and WriteStream() functions.  Streaming of
  base classes must be done by explicitly calling base::ReadStream() or base::WriteStream().
  There is no special support for streaming virtual base classes.  

  Versioning support:

  This is designed to allow versioning of sub-objects, without
  each sub-object needing to keep track of version numbers itself.
  For example, suppose we have a class foo that has two versions,
  version 1 and version 2.  Both versions contain a sub-object bar,
  which has a different on-disk format between the two versions,
  which doesn't allow for a separate version number.  Then
  the serializer for foo will register with the pstream object
  a version number that the bar serializer will be able to read.

  Version numbers are tagged by a value of type VersionTag,
  and are stored in a stack structure inside the pstream object.
  For example, to indicate that sub-objects should expect version 2 for
  the tag FunctionDatabaseVersion, a serializer will be defined as:

  // global tag object, with default version number 1
  PStream::VersionTag FunctionDatabaseV(1);

  ipstream& operator>>(ipstream& in, foo& f)
  {
     PStream::VersionSentry Sentry(in, FunctionDatabaseVersion , 2);
     in >> f.b;
     return in;
  }

  If the version number is to be read from the stream, then 
  use in.read<int>() as the version number.

  While the VersionTag object is in scope, then the pstream will keep
  track of the associated version number, which can be read from the pstream object:

  int Version = in.version_of(FunctionDatabaseV);

  The VersionSentry also tracks the version number, as Sentry.version()

  Alternatively, the version number can be set manually, by
  in.push_version(FunctionDatabaseV, 2);
  and reset with in.pop_version(FunctionDatabaseV);
  But this method is not recommended, except to set initial version numbers
  on a new stream, if they differ from the default.

  Experimental features not yet merged: streaming of typemap's.
*/

#if !defined(PSTREAM_H_437Y8RFH838Y89FWYQ378O4YH347YW78)
#define PSTREAM_H_437Y8RFH838Y89FWYQ378O4YH347YW78

#include "pstreamfwd.h"
#include "common/trace.h"
#include "common/ctassert.h"
#include <boost/type_traits.hpp>
#include <iostream>
#include <utility>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <stack>
#include <complex>

namespace PStream
{

// see pstreamformats.h for {i|o}pstreambuf formats

int const NumConcreteStreamClasses = 4;

template <int Format>
class opstreambuf;

template <int Format>
class ipstreambuf;

class generic_opstreambuf;
class generic_ipstreambuf;

template <int Format>
struct pstreambuf_traits;

template <typename T>
struct pstream_type_traits
{
   static bool const is_fundamental = boost::is_fundamental<T>::value;
};

class VersionTag;

//
// opstream
//

class opstream
{
   public:
      virtual ~opstream();

      // returns the unique ID of the format
      int format() const { return BaseFormat; }

      // writes an id to the stream.  Default version raises an exception/abort
      virtual void put_id(id_type Id);

      virtual void flush();  // default version calls overflow()

      generic_opstreambuf* get_buffer() const { return Buffer; }

      template <typename T>
      void write(T const& x);

      void push_version(VersionTag const& Tag, int Version);
      void pop_version(VersionTag const& Tag);
      int version_of(VersionTag const& Tag) const;

   protected:
      opstream(int Format_, generic_opstreambuf* Buffer_);

      opstream(int Format_, byte_type* BufBegin, byte_type* BufEnd, byte_type* BufPtr);

      virtual void overflow() = 0;

      // changes the format of the stream to OtherFormat
      void set_format(int OtherFormat);

      // changes the format of the stream to OtherFormat and writes the format ID
      // to the stream.
      void put_format(int OtherFormat);

      // writes the format ID to the stream
      void put_format();

      void set_buffer(byte_type* BufBegin, byte_type* BufEnd, byte_type* BufPtr);

      byte_type* buf_begin() const;
      byte_type* buf_end() const;
      byte_type* buf_ptr() const;

      void set_buf_ptr(byte_type* Ptr);

   private:
      opstream(opstream&);                  // not implemented
      opstream& operator=(opstream const&); // not implemented

      int BaseFormat;
      generic_opstreambuf* Buffer;

      typedef std::stack<int> VersionStackType;
      std::map<VersionTag const*, VersionStackType> VersionNumbers;

   friend class generic_opstreambuf; // so that generic_opstream can call overflow()
};

class generic_opstreambuf
{
   public:
      virtual ~generic_opstreambuf();

      void put(byte_type b);

      void put_n(byte_type const* Src, size_t Count);

      template <int Atom>
      void put_bswap(byte_type const* Src);

      template <int Atom>
      void put_n_bswap(byte_type const* Src, size_t Count);

      //      void put_n_reverse(byte_type const* Src, size_t Count);
      
      // direct access to the buffer
      byte_type* buf_begin() const { return BufBegin; }
      byte_type* buf_end() const { return BufEnd; }
      byte_type* buf_ptr() const { return BufPtr; }

      void set_buf_ptr(byte_type* Ptr) { BufPtr = Ptr; }

      // returns the number of bytes available
      size_t buf_avail() const { return BufEnd - BufPtr; }

      // increments the buf_ptr by the specified count, for direct writing to the buffer.
      void bump(size_t Count);

      void flush();

      void set_buffer(byte_type* Begin, byte_type* End, byte_type* Ptr);

      operator opstream&() { return *BackPointer; }

      void put_id(id_type Id) { BackPointer->put_id(Id); }

      void push_version(VersionTag const& Tag, int Version);
      void pop_version(VersionTag const& Tag);
      int version_of(VersionTag const& Tag);

   protected:
      generic_opstreambuf(opstream* BackPointer_,
			  byte_type* Beg, byte_type* End, byte_type* Ptr);

      generic_opstreambuf(generic_opstreambuf const& Buffer);

      void overflow();

      opstream* BackPointer;

      byte_type* BufBegin;
      byte_type* BufEnd;
      byte_type* BufPtr;

   private:
      static generic_opstreambuf* CopyOpstreambuf(int Format, generic_opstreambuf const& Buffer);
      static generic_opstreambuf* MakeOpstreambuf(int Format, 
						  opstream* BackPointer,
						  byte_type* BufBegin, 
						  byte_type* BufEnd, 
						  byte_type* BufPtr);

      template <int Atom>
      void put_bswap_unchecked(byte_type const* Src);

      generic_opstreambuf& operator=(generic_opstreambuf const&); // not implemented

   friend class opstream;
};

template <int Format>
class opstreambuf : public generic_opstreambuf
{
   public:
      typedef typename pstreambuf_traits<Format>::size_type       size_type;
      typedef typename pstreambuf_traits<Format>::difference_type difference_type;

      template <typename T>
      void write(T const& x);

      template <typename T>
      void write_n(T const* x, size_t N);

   private:
      ~opstreambuf();

      opstreambuf(opstream* BackPointer_,
		  byte_type* Beg, byte_type* End, byte_type* Ptr);


      opstreambuf(generic_opstreambuf const& Other)
	: generic_opstreambuf(Other) {}

      opstreambuf(opstreambuf&);             // not implemented
      opstreambuf& operator==(opstreambuf&); // not implemented

   friend class generic_opstreambuf;
};

//
// ipstream
//

class ipstream
{
   public:
      virtual ~ipstream();

      // returns the unique ID of the format
      int format() const { return BaseFormat; }

      // reads an id from the stream.  Default version raises an exception/abort
      virtual id_type get_id();

      generic_ipstreambuf* get_buffer() const { return Buffer; }

      template <typename T>
      T read();

      void push_version(VersionTag const& Tag, int Version);
      void pop_version(VersionTag const& Tag);
      int version_of(VersionTag const& Tag) const;

   protected:
      ipstream(int Format_, generic_ipstreambuf* Buffer_);

      ipstream(int Format_, byte_type const* BufBegin, byte_type const* BufEnd, byte_type const* BufPtr);

      virtual void underflow() = 0;

      // changes the format of the stream to OtherFormat
      void set_format(int OtherFormat);

      // reads the format ID from the stream and changes to it.  Returns the new format.
      int get_format();

      void set_buffer(byte_type const* BufBegin, 
		      byte_type const* BufEnd, 
		      byte_type const* BufPtr);

      byte_type const* buf_begin() const;
      byte_type const* buf_end() const;
      byte_type const* buf_ptr() const;

      void set_buf_ptr(byte_type const* Ptr);
      void set_buf_end(byte_type const* Ptr);

   private:
      ipstream(ipstream&);                  // not implemented
      ipstream& operator=(ipstream const&); // not implemented

      int BaseFormat;
      generic_ipstreambuf* Buffer;

      typedef std::stack<int> VersionStackType;
      std::map<VersionTag const*, VersionStackType> VersionNumbers;

   friend class generic_ipstreambuf;
};

class generic_ipstreambuf
{
   public:
      byte_type get();

      void get_n(byte_type* Dest, size_t Count);

      void get_n_reverse(byte_type* Dest, size_t Count);

      // direct access to the buffer
      byte_type const* buf_begin() const { return BufBegin; }
      byte_type const* buf_end() const { return BufEnd; }
      byte_type const* buf_ptr() const { return BufPtr; }

      // returns the number of bytes available before underflow
      size_t buf_avail() const { return BufEnd - BufPtr; }

      // increments the buf_ptr by the specified count, for direct reading from the buffer.
      void bump(size_t Count);

      void set_buf_ptr(byte_type const* Ptr) { BufPtr = Ptr; }
      void set_buf_end(byte_type const* End) { BufEnd = End; }

      void set_buffer(byte_type const* Begin, byte_type const* End, byte_type const* Ptr);

      operator ipstream&() { return *BackPointer; }

      id_type get_id() { return BackPointer->get_id(); }

      void push_version(VersionTag const& Tag, int Version);
      void pop_version(VersionTag const& Tag);
      int version_of(VersionTag const& Tag);

   protected:
      generic_ipstreambuf(ipstream* BackPointer_,
			  byte_type const* Beg, byte_type const* End, byte_type const* Ptr);

      virtual ~generic_ipstreambuf();

      generic_ipstreambuf(generic_ipstreambuf const& Buffer);

      void underflow();

      ipstream* BackPointer;

      byte_type const* BufBegin;
      byte_type const* BufEnd;
      byte_type const* BufPtr;

   private:
      static generic_ipstreambuf* CopyIpstreambuf(int Format, generic_ipstreambuf const& Buffer);
      static generic_ipstreambuf* MakeIpstreambuf(int Format, 
						  ipstream* BackPointer,
						  byte_type const* BufBegin, 
						  byte_type const* BufEnd, 
						  byte_type const* BufPtr);

      generic_ipstreambuf& operator=(generic_ipstreambuf const&); // not implemented

   friend class ipstream;
};

template <int Format>
class ipstreambuf : public generic_ipstreambuf
{
   public:
      typedef typename pstreambuf_traits<Format>::size_type       size_type;
      typedef typename pstreambuf_traits<Format>::difference_type difference_type;

      template <typename T>
      T read();

   private:
      ~ipstreambuf();

      ipstreambuf(ipstream* BackPointer_,
		  byte_type const* Beg, byte_type const* End, byte_type const* Ptr);


      ipstreambuf(generic_ipstreambuf const& Other)
	: generic_ipstreambuf(Other) {}

      ipstreambuf(ipstreambuf&);             // not implemented
      ipstreambuf& operator==(ipstreambuf&); // not implemented

   friend class ipstream;
   friend class generic_ipstreambuf;
};

//
// opstream_iterator
//

template <typename T>
class opstream_iterator : public std::iterator<std::output_iterator_tag, T, void, void, void>
{
   public:
      opstream_iterator(opstream& out_) : out(&out_) {}
      opstream_iterator(opstream_iterator const& os) : out(os.out) {}

      opstream_iterator& operator=(T const& value)
      {
 	 *out << value;
	 return *this;
      }

      opstream_iterator& operator*()      { return *this; }
      opstream_iterator& operator++()     { return *this; } 
      opstream_iterator& operator++(int)  { return *this; } 

 private:
      opstream* out;
};

template <int Format, typename T>
class opstreambuf_iterator : public std::iterator<std::output_iterator_tag, T, void, void, void>
{
   public:
      opstreambuf_iterator(opstreambuf<Format>& out_) : out(&out_) {}
      opstreambuf_iterator(opstreambuf_iterator<Format, T> const& os) : out(os.out) {}

      opstreambuf_iterator& operator=(T const& value) 
      {
	 *out << value;
	 return *this;
      }

      opstreambuf_iterator& operator*()      { return *this; }
      opstreambuf_iterator& operator++()     { return *this; } 
      opstreambuf_iterator& operator++(int)  { return *this; } 

 private:
      opstreambuf<Format>* out;
};

/*
  We could have an overload of std::copy() (although we cannot put it in namespace std)
  that unwraps copying to an opstream_iterator to a copy to an opstreambuf_iterator<Format>.
  However we don't know if this will be efficient; if there is no opstreambuf<Format> inserter
  for the type and it does the conversion back to opstream for every operation then it will
  be slightly inefficient.  Instead, we leave the choice up to the caller.  If you want
  an iterator operation to be efficient, then make sure you are using an opstreambuf_iterator.
*/

//
// ipstream_iterator
//

template <typename T>
class ipstream_iterator : public std::iterator<std::input_iterator_tag, T, ptrdiff_t, T const*, T const&>
{
   public:
      ipstream_iterator(ipstream& in_) : in(&in_) {}
      ipstream_iterator(ipstream_iterator const& is) : in(is.in) {}

      T const& operator*() const { (*in) >> value; return value; }
      T const* operator->() const { (*in) >> value; return &value; }

      ipstream_iterator& operator++()     { return *this; }
      ipstream_iterator& operator++(int)  { return *this; }

 private:
      ipstream* in;
      mutable T value;
};

template <int Format, typename T>
class ipstreambuf_iterator : public std::iterator<std::input_iterator_tag, T, ptrdiff_t, T const*, T const&>
{
   public:
      ipstreambuf_iterator(ipstreambuf<Format>& in_) : in(&in_) {}
      ipstreambuf_iterator(ipstreambuf_iterator const& is) : in(is.in) {}

      T const& operator*() const { (*in) >> value; return value; }
      T const* operator->() const { (*in) >> value; return &value; }

      ipstreambuf_iterator& operator++()     { return *this; } 
      ipstreambuf_iterator& operator++(int)  { return *this; } 

 private:
      ipstreambuf<Format>* in;
      mutable T value;
};

// Versioning support

class VersionTag
{
   public:
      VersionTag(); // not implemented
      VersionTag(VersionTag const&); // not implemented
      VersionTag& operator=(VersionTag const&); // not implemented

      VersionTag(int DefaultVersion_) : DefaultVersion(DefaultVersion_) {}

      int default_version() const { return DefaultVersion; }

   private:
      int DefaultVersion;
};

class VersionSentry
{
   public:
      VersionSentry(); // not implemented
      VersionSentry(VersionSentry const&); // not implemented
      VersionSentry& operator=(VersionSentry const&); // not implemented

      VersionSentry(ipstream& In, VersionTag const& Tag_, int v_);

      template <int Format>
      VersionSentry(ipstreambuf<Format>& In, VersionTag const& Tag_, int v_);

      VersionSentry(opstream& Out, VersionTag const& Tag_, int v_);

      template <int Format>
      VersionSentry(opstreambuf<Format>& Out, VersionTag const& Tag_, int v_);

      ~VersionSentry();

      int version() const { return v; }
      
   private:
      VersionTag const& Tag;
      ipstream* in;
      opstream* out;
      int v;
};

//
// the main inserter/extractor.
//
// These only get called if there is no better overload.
// This implements the single-dispatch on Buf.which()
// delegates to operator<<() on the derived stream type.
//

template <typename T>
opstream& operator<<(opstream& Buf, T const& Obj);

template <typename T>
ipstream& operator>>(ipstream& Buf, T& Obj);

//
// Support for streaming polymorphic objects.
// The underlying library is capable of streaming polymorphic objects
// from an arbitary heirachy (no particular need to inherit from StreamableBase),
// however extra support is required, specifically a generic version
// of StreamExtract().  Maybe this is actually possible?  Does ConstructFromStream work for this?
// In any event, it is still convenient to have a StreamableBase class.
//

//
// poly_traits
//
// Determines whether a given type is polymorphic or not.  This is probably obsoleted
// by boost, and almost certainly will be obsoleted by C++0x.
// This implementation uses the non-portable hack of checking whether an additional
// virtual method changes the size of the type.  If the type is not polymorphic,
// then a virtual function introduces a vtable and sizeof() grows by a pointer.
// if the type is alreadypolymorphic, then the sizeof() does not change.
// There are probably C++ compilers out there where this does not work.
// A special override is necessary for builtins because it is not possible to derive from them.
//

template <typename T>
struct poly_traits
{
   private:
      class TestPoly : public T
      {
         public:
	    virtual ~TestPoly();
      };

   public:
      static const bool is_polymorphic = (sizeof(T) == sizeof(TestPoly));
};

#define DEFINE_BUILTIN_POLY_TRAITS(type)	\
template <>					\
struct poly_traits<type>			\
{						\
   static const bool is_polymorphic = false;	\
};

DEFINE_BUILTIN_POLY_TRAITS(bool)
DEFINE_BUILTIN_POLY_TRAITS(char)
DEFINE_BUILTIN_POLY_TRAITS(unsigned char)
DEFINE_BUILTIN_POLY_TRAITS(signed char)
DEFINE_BUILTIN_POLY_TRAITS(short)
DEFINE_BUILTIN_POLY_TRAITS(unsigned short)
DEFINE_BUILTIN_POLY_TRAITS(int)
DEFINE_BUILTIN_POLY_TRAITS(unsigned int)
DEFINE_BUILTIN_POLY_TRAITS(long)
DEFINE_BUILTIN_POLY_TRAITS(unsigned long)
#if defined(USE_LONGLONG)
DEFINE_BUILTIN_POLY_TRAITS(long long)
DEFINE_BUILTIN_POLY_TRAITS(unsigned long long)
#endif
DEFINE_BUILTIN_POLY_TRAITS(float)
DEFINE_BUILTIN_POLY_TRAITS(double)

template <typename T>
struct poly_traits<T*>
{
   static bool const is_polymorphic = false;
};

//
// streamable_traits
//
// default version uses the non-portable hack to determine is_polymorphic
// and assumes clonable semantics.  Specialize streamable_traits for
// your type if the defaults are not correct.
//

// generic case covers is_poly == true
template <typename T, bool is_poly>
struct streamable_traits_helper
{
   static bool const is_polymorphic = true;
   static bool const is_clonable = true;
   static T* clone(T const* Obj) { return Obj->clone(); }
};

// partial specialization for is_poly == false
template <typename T>
struct streamable_traits_helper<T, false>
{
   static bool const is_polymorphic = false;
   static bool const is_clonable = true;
   static T* clone(T const* Obj) { return new T(Obj); }
};

template <typename T>
struct streamable_traits : public streamable_traits_helper<T, poly_traits<T>::is_polymorphic>
{
};

class StreamableBase
{
   public:
      StreamableBase() {}
      virtual ~StreamableBase() {}

      virtual void WriteStream(opstream& stream) const = 0;
      virtual void ReadStream(ipstream& in); // the default implementation does nothing

      static StreamableBase* ConstructFromStream(ipstream& in); // "default" returns NULL

      virtual StreamableBase* clone() const { return NULL; }
};

class StreamableFactory
{
   public:
      StreamableFactory() {}
      virtual ~StreamableFactory() {}

      virtual StreamableBase* construct(ipstream& stream) = 0;
};

// Registers a class for persistent streaming, with the unique identifier 's'
template <typename T>
void Register(std::string const& s);

StreamableBase* StreamExtract(ipstream& stream);

void RegisterClass(const std::type_info& t, const std::string& UID, StreamableFactory* Factory);

// GetUID returns the UID of the given type.  If the type is not registered,
// then GetUID throws an exception
std::string GetUID(const std::type_info& t);

// CheckForUID returns the UID of the given type.  If the type is not registered,
// then CheckForUID returns the empty string.
std::string CheckForUID(const std::type_info& t);

StreamableBase* StreamExtract(ipstream& stream);

template <typename T, bool IsPolymorphic = poly_traits<T>::is_polymorphic>
struct ConstructHelper;

template <typename T>
struct ConstructHelper<T, false>
{
   static T* apply(ipstream& in)
   {
      T* Result = new T;
      in >> *Result;
      return Result;
   }
};

template <typename T>
struct ConstructHelper<T, true>
{
   static T* apply(ipstream& in)
   {
      return &dynamic_cast<T&>(*StreamExtract(in));
   }
};

template <typename T>
T* ConstructFromStream(ipstream& in)
{
   return ConstructHelper<T>::apply(in);
}

//
// streaming of STL types
//

template <int Format>
opstreambuf<Format>& operator<<(opstreambuf<Format>& out, char const* String);

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::basic_string<T> const& vec);

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::basic_string<T>& vec);

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::vector<T> const& vec);

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::vector<T>& vec);

// std::vector<bool> is not a standard container, and requires special case code
// to stream it.
template <int Format>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::vector<bool> const& vec);

template <int Format>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::vector<bool>& vec);

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::deque<T>& vec);

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::deque<T>& vec);

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::list<T>& vec);

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::list<T>& vec);

template <int Format, typename T, class Comp>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::set<T, Comp> const& vec);

template <int Format, typename T, class Comp>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::set<T, Comp>& vec);

template <int Format, typename T1, typename T2, class Comp>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, const std::map<T1, T2, Comp>& vec);

template <int Format, typename T1, typename T2, class Comp>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::map<T1, T2, Comp>& vec);

template <int Format, typename T1, typename T2>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::pair<T1, T2> const& Obj);

template <int Format, typename T1, typename T2>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::pair<T1, T2>& Obj);

template <int Format, typename T>
opstreambuf<Format>& operator<<(opstreambuf<Format>& stream, std::complex<T> const& vec);

template <int Format, typename T>
ipstreambuf<Format>& operator>>(ipstreambuf<Format>& stream, std::complex<T>& vec);

} // namespace pstream

#include "pstream.cc"

#endif
