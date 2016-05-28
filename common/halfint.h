// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/halfint.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(HALFINT_H_FDS78347FVUFUIJ89UEJI389URVJI)
#define HALFINT_H_FDS78347FVUFUIJ89UEJI389URVJI

#include <iostream>
#include <stdexcept>
#include <math.h>
#include <string>
#include "convertstring.h"
#include "trace.h"

#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

class half_int
{
   public:
      half_int() : N2(0) {}  // initialize to zero, rather than leave N2 unitialized.
      half_int(double D) { N2 = int(floor(D * 2.0 + 0.5)); }     // round 2D to nearest int
      half_int(int I) { N2 = 2 * I; }

      // use compiler defined copy ctor, assignment and dtor
      // (user defined ctor affects performance on intel 6 compiler)

      struct twice_tag {};
      half_int(int N, twice_tag) : N2(N) {}

      const half_int& operator=(const half_int& H) { N2 = H.N2; return *this; }

      half_int& operator+=(const half_int& H) { N2 += H.N2; return *this; }
      half_int& operator-=(const half_int& H) { N2 -= H.N2; return *this; }

      half_int& operator++() { N2 += 2; return *this; }
      half_int& operator--() { N2 -= 2; return *this; }

      half_int operator++(int) { half_int Temp(*this); N2 += 2; return Temp; }
      half_int operator--(int) { half_int Temp(*this); N2 -= 2; return Temp; }

      double to_double() const { return N2 / 2.0; }

      int to_int() const 
	{ if (!is_integral()) throw_cannot_convert(); return N2 >> 1; }

      // fast version of to_int, where is_integral() is a precondition
      int to_int_assert() const 
	{ DEBUG_PRECONDITION(is_integral()); return N2 >> 1; }

      bool is_integral() const { return !(N2 & 1); }

      int twice() const { return N2; }

   private:
      int N2;

      // inline throws are expensive on KAI, so delegate to a non-inlined function
      static void throw_cannot_convert();

#if defined(USE_PSTREAM)
   friend PStream::opstream& operator<<(PStream::opstream& out, const half_int& H);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, half_int& H);
#endif

};

inline
half_int from_twice(int t)
{
   return half_int(t, half_int::twice_tag());
}

// non-members acting on a half-int

bool operator==(half_int h1, half_int h2);
bool operator!=(half_int h1, half_int h2);

bool operator<(half_int h1, half_int h2);
bool operator>(half_int h1, half_int h2);
bool operator<=(half_int h1, half_int h2);
bool operator>=(half_int h1, half_int h2);

half_int operator*(half_int h, int i);
half_int operator*(int i, half_int h);
double   operator*(half_int h, double d);
double   operator*(double d, half_int h);
double   operator*(half_int h1, half_int h2);

double   operator/(half_int h, double d);

half_int operator+(half_int h1, half_int h2);
half_int operator-(half_int h1, half_int h2);

half_int abs(half_int h);

half_int operator-(half_int h);

bool is_integral(half_int h);

int to_int(half_int h);             // throws if h is not integral
int to_int_assert(half_int h);       // assumes h is integral as a precondition (faster)

// conversion to string, as a fraction,
// if h is half-integral then the string is of the forn "n/2"
std::string to_string_fraction(half_int h);

// A function for (-1)^x
long minus1pow(long x);
int minus1pow(int x);

std::ostream& operator<<(std::ostream& out, const half_int& H);
std::istream& operator>>(std::istream& in, half_int& H);

// convert_string specialization
// Currently, this doesn't detect overflow in some circumstances.  It would only affect
// half_int's greater than numeric_limits<int>::max() / 2, which presumably never happen ;)
template <class FwdIter>
struct convert_string_partial<half_int, FwdIter>
{
   static half_int apply(FwdIter start, FwdIter end);
};


// inlines

inline
bool operator==(half_int h1, half_int h2)
{
   return h1.twice() == h2.twice();
}

inline
bool operator!=(half_int h1, half_int h2)
{
   return h1.twice() != h2.twice();
}

inline
bool operator<(half_int h1, half_int h2)
{
   return h1.twice() < h2.twice();
}

inline
bool operator>(half_int h1, half_int h2)
{
   return h1.twice() > h2.twice();
}

inline
bool operator<=(half_int h1, half_int h2)
{
   return h1.twice() <= h2.twice();
}

inline
bool operator>=(half_int h1, half_int h2)
{
   return h1.twice() >= h2.twice();
}

inline
half_int operator*(half_int h, int i)
{
   return half_int(h.twice() * i, half_int::twice_tag());
}

inline
half_int operator*(int i, half_int h)
{
   return half_int(i * h.twice(), half_int::twice_tag());
}

inline
double operator*(half_int h, double d)
{
   return h.twice() * d / 2.0;
}

inline
double operator*(double d, half_int h)
{
   return d * h.twice() / 2.0;
}

inline
double operator*(half_int h1, half_int h2)
{
   return h1.twice() * h2.twice() / 4.0;
}

inline
double operator/(half_int h, double d)
{
   return h.twice() / (d * 2.0);
}

inline
half_int operator+(half_int h1, half_int h2)
{
   return half_int(h1.twice() + h2.twice(), half_int::twice_tag());
}

inline
half_int operator-(half_int h1, half_int h2)
{
   return half_int(h1.twice() - h2.twice(), half_int::twice_tag());
}

inline
half_int abs(half_int h)
{
   return half_int(abs(h.twice()), half_int::twice_tag());
}

inline
half_int operator-(half_int h)
{
   return half_int(-h.twice(), half_int::twice_tag());
}

inline 
bool is_integral(half_int h)
{
   return h.is_integral();
}

inline 
int to_int(half_int h)
{
   return h.to_int();
}

inline 
int to_int_assert(half_int h)
{
   return h.to_int_assert();
}

// there isn't a more sensible place to define these?
inline
long minus1pow(long x)
{
   return 1 - ((x & 1) << 1);
}

inline
int minus1pow(int x)
{
   return 1 - ((x & 1) << 1);
}

//
// returns true if (a,b,c) satisfies the triangle condition 
// |a-c| <= b <= a+c and a+b+c is integral.
// a,b,c are assumed positive.
inline
bool is_triangle(half_int a, half_int b, half_int c)
{
   return (b >= abs(a-c)) && (b <= a+c) && (a+b+c).is_integral();
}

#if defined(USE_PSTREAM)
inline
PStream::opstream& operator<<(PStream::opstream& out, const half_int& H)
{
   return out << H.N2;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, half_int& H)
{
   return in >> H.N2;
}
#endif // #defined (USE_PSTREAM)

#include "halfint.cc"

#endif
