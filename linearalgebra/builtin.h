/* -*- C++ -*- $Id$
  Traits to return the resulting type of a "usual arithmetic conversion";
  see section 5.2 of the C++ standard.

  Created 2004-04-09 Ian McCulloch
*/

#if !defined(BUILTIN_H_SADHCIUEWHRIU43Y587Y87HYDO8HW874HYOHFWAO)
#define BUILTIN_H_SADHCIUEWHRIU43Y587Y87HYDO8HW874HYOHFWAO

namespace LinearAlgebra
{

//
// builtin_conversion
//
// From the C++ Standard, section 5.2:
//    Many binary operators that expect operands of arithmetic type cause 
//    conversions and yield result types in a similar way. The purpose is 
//    to yield a common type, which is also the type of the result. This 
//    pattern is called the usual arithmetic conversions. 
//
// For builtin arithmetic types, builtin_conversion<T1, T2>::type
// is the type of the expression resulting from "usual arithmetic conversion",
// defined in section 5.2 of the C++ standard.
//

template <typename T1, typename T2>
struct builtin_conversion;

namespace Private
{

template <typename T> struct type_weight;

// conversions of wchar_t are not included here.

template <> struct type_weight<bool> { static int const value = 0; };

// This will fail if char is unsigned && sizeof(char) == sizeof(int).
template <> struct type_weight<char> { static int const value = 0; };

template <> struct type_weight<signed char> { static int const value = 0; };
template <> struct type_weight<unsigned char> 
{ static int const value = sizeof(unsigned char) >= sizeof(int); };
template <> struct type_weight<short> { static int const value = 0; };
template <> struct type_weight<unsigned short> 
{ static int const value = sizeof(unsigned short) >= sizeof(int); };
template <> struct type_weight<int> { static int const value = 0; };

template <> struct type_weight<unsigned int> { static int const value = 1; };
template <> struct type_weight<long> { static int const value = 2; };
template <> struct type_weight<unsigned long> { static int const value = 3; };
#if defined(USE_LONGLONG)
template <> struct type_weight<long long> { static int const value = 4; };
template <> struct type_weight<unsigned long long> { static int const value = 5; };
#endif
template <> struct type_weight<float> { static int const value = 7; };
template <> struct type_weight<double> { static int const value = 9; };

/*
  If the highest weight type is an unsigned type, or if
  both types are signed, then the result type is equal to the highest weight type.

  The last case is tricky: the hightest weight type is signed and the
  lowest weight type is unsigned.  If the heighest weight (signed) type can
  represent all values of the lower type (which I interpret to mean
  strict inequality sizeof(lowest) < sizeof(highest)), then the result
  type is the hightest weight (signed) type.  Otherwise, the highest
  weight type is 'promoted' to an unsigned type.
*/

template <typename T1, typename T2, 
	  int UnsignedSigned = (type_weight<T1>::value % 2 != 0) &&
	  (type_weight<T2>::value % 2 == 0)>
struct promote_unsigned;

template <typename T1, typename T2>
struct promote_unsigned<T1, T2, 0>
{
   static int const value = type_weight<T2>::value;
};

template <typename T1, typename T2>
struct promote_unsigned<T1, T2, 1>
{
   static int const value = type_weight<T2>::value + (sizeof(T1) == sizeof(T2));
};

template <int weight>
struct type_from_weight { static int const value = -1; };

template <> struct type_from_weight<0>  { typedef int type; };
template <> struct type_from_weight<1>  { typedef unsigned int type; };
template <> struct type_from_weight<2>  { typedef long type; };
template <> struct type_from_weight<3>  { typedef unsigned long type; };
#if defined(USE_LONGLONG)
template <> struct type_from_weight<4>  { typedef long long type; };
template <> struct type_from_weight<5>  { typedef unsigned long long type; };
#endif
template <> struct type_from_weight<7>  { typedef float type; };
template <> struct type_from_weight<9>  { typedef double type; };

} // namespace Private

template <typename T1, typename T2>
struct builtin_conversion
{
   private:
      static int const v1 = Private::type_weight<T1>::value;
      static int const v2 = Private::type_weight<T2>::value;
      static int const value = v1 <= v2 ? Private::promote_unsigned<T1, T2>::value
				   : Private::promote_unsigned<T2, T1>::value;
   public:
      typedef typename Private::type_from_weight<value>::type type;
};

} // namespace LinearAlgebra

#endif
