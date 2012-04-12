/* -*- C++ -*- $Id$

  Defines appropriate typedef's for the possible binary formats
  used by PStreams.

  Created 2002-04-08 Ian McCulloch


  Synopsis:

  #include "endian.h"
  namespace CurrentFormat
  {
     const EndianType Endianness = <implementation defined>;

     typedef <implementation defined> char8;
     typedef <implementation defined> uchar8;
     typedef <implementation defined> int16;
     typedef <implementation defined> uint16;
     typedef <implementation defined> int32;
     typedef <implementation defined> uint32;
     typedef <implementation defined> int64;
     typedef <implementation defined> uint64;
 
     typedef <implementation defined> ieee32;
     typedef <implementation defined> ieee64;

     template <typename Builtin, typename FormatType> 
     struct IsTrivialConversion
     {
        static bool const value = <implementation defined>;
     };
  }

  Requirements:

  The namespace CurrentFormat defines typedef's mapping
  the chosen binary types onto builtin (or possibly user-defined) types.  
  Currently, we assume that the representation and size of the builtin type
  exactly corresponds with the external representation of the given endianness.
  On some platforms that we almost certainly don't care about
  this may not be possible to achieve; in that case an extra indirection would
  be required to do the appropriate conversion,
  somewhere between CurrentFormat and the actual stream I/O functions.

  The list of binary types is:

  Integral types (signed types are 2's compliment):

  char8           signed    8 bit                   -128   ..                   127
  uchar8          unsigned  8 bit                      0   ..                   255
  int16           signed   16 bit                 -32768   ..                 32767
  uint16          unsigned 16 bit                      0   ..                 65536
  int32           signed   32 bit            -2147483648   ..            2147483647
  uint32          unsigned 32 bit                      0   ..            4294967295
  int64           signed   64 bit   -9223372036854775808   ..   9223372036854775807
  uint64          unsigned 64 bit                      0   ..  18446744073709551615

  Floating point types (IEEE-754-1985)

  ieee32          IEEE single-precision
  ieee64          IEEE double-precision

  The platform must have a definite endianness (ie, big endian or little endian),
  which is defined in CurrentFormat::Endianness.

  The platform must also define a traits type 
  CurrentFormat::IsTrivialConversion<typename Builtin, typename FormatType>,
  where Intrinsic is a builtin type, FormatType is one of the above typedefs.
  This traits type must contain a boolean constant 'value' which is set to
  true if a given builtin type has a trivial conversion to the FormatType 
  (meaning a memcpy() or reinterpret_cast will work, with no possibility 
  of overflow or underflow).

  IsTrivialConversion doesn't test endianness; it is used to determine only if
  the given builtin type can be converted trivially to the representation of
  a given format type using the host endianness.  

  If IsTrivialConversion is true, and additionally
  the endianness of the target format is the same as the host endianness, then
  a memcpy() could be used to stream values of that builtin type.  Otherwise,
  if IsTrivialConversion is true but the endianness doesn't match, then
  an endianness conversion (but no representation conversion) will be required.
  If IsTrivialConversion is false, then a conversion of representation is required,
  perhaps using boost::numeric_cast<>.

  Builtin type 'bool' maps onto int8.
*/

#if !defined(CURRENTFORMAT_H_YEUWIHFE67834Y723OY723OY23824398423Y)
#define CURRENTFORMAT_H_YEUWIHFE67834Y723OY723OY23824398423Y

#include "config.h"
#include "endian.h"
#include <boost/cstdint.hpp>
namespace CurrentFormat
{

#if defined(WORDS_BIGENDIAN)
const EndianType Endianness = BigEndian;
#else
const EndianType Endianness = LittleEndian;
#endif // defined(WORDS_BIGENDIAN)

typedef boost::int8_t char8;
typedef boost::uint8_t uchar8;

typedef boost::int16_t int16;
typedef boost::uint16_t uint16;

typedef boost::int32_t int32;
typedef boost::uint32_t uint32;

typedef boost::int64_t int64;
typedef boost::uint64_t uint64;

#if !defined(SIZEOF_FLOAT) || !defined(SIZEOF_DOUBLE) || !defined(SIZEOF_INT) || !defined(SIZEOF_LONG)
#error "SIZEOF_XXX macros are not defined in config.h, problem with configure script?"
#endif

#if SIZEOF_FLOAT != 4
#error "float is not a 32-bit IEEE type!"
#endif
typedef float ieee32;

#if SIZEOF_DOUBLE != 8
#error "double is not a 64-bit IEEE type!"
#endif
typedef double ieee64;

template <typename Builtin, typename FormatType>
struct IsTrivialConversion
{
   static bool const value = false;
};

struct Trivial { static bool const value = true; };

template <typename T> struct IsTrivialConversion<T, T> : public Trivial {};

template <> struct IsTrivialConversion<bool, char8> : public Trivial {};

template <> struct IsTrivialConversion<signed char, char8> : public Trivial {};
template <> struct IsTrivialConversion<unsigned char, uchar8> : public Trivial {};
#if defined(__CHAR_UNSIGNED__)
template <> struct IsTrivialConversion<char, uchar8> : public Trivial {};
#else
template <> struct IsTrivialConversion<char, char8> : public Trivial {};
#endif

#if SIZEOF_SHORT == 2
template <> struct IsTrivialConversion<short, int16> : public Trivial {};
template <> struct IsTrivialConversion<unsigned short, uint16> : public Trivial {};
#endif

#if SIZEOF_INT == 4
template <> struct IsTrivialConversion<int, int32> : public Trivial {};
template <> struct IsTrivialConversion<unsigned int, uint32> : public Trivial {};
#endif

#if SIZEOF_LONG == 4
template <> struct IsTrivialConversion<long, int32> : public Trivial {};
template <> struct IsTrivialConversion<unsigned long, uint32> : public Trivial {};
#elif SIZEOF_LONG == 8
template <> struct IsTrivialConversion<long, int64> : public Trivial {};
template <> struct IsTrivialConversion<unsigned long, uint64> : public Trivial {};
#endif

#if defined(USE_LONGLONG)
# if SIZEOF_LONG_LONG == 8
template <> struct IsTrivialConversion<long long, int64> : public Trivial {};
template <> struct IsTrivialConversion<unsigned long long, uint64> : public Trivial {};
# endif
#endif

} // namespace CurentFormat

#endif
