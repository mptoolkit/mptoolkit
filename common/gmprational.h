// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/gmprational.h
//
// Copyright (C) 2001-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  gmprational.h

  Wraps the GNU-MP library for rationals.

  Created 2001-04-01
*/

#if !defined(GMPRATIONAL_H_SFDHU437Y8EH7REHUG987YFHURFIU4738YR4)
#define GMPRATIONAL_H_SFDHU437Y8EH7REHUG987YFHURFIU4738YR4

#include "gmp.h"
#include "gmpint.h"

namespace gmp
{

class rational
{
   public:
      rational();
      rational(rational const& x);
      rational(int x);
      rational(long x);
      rational(unsigned long x);
      rational(bigint const& x);
      rational(double x);        // exact, no rounding

      explicit rational(char const* x);            // not yet implemented
      explicit rational(std::string const& x);     // not yet implemented

      ~rational();

      rational& operator=(rational const& x);
      rational& operator=(int x);
      rational& operator=(long x);
      rational& operator=(unsigned long x);

      void swap(rational& x);

      rational& operator+=(rational const& x);
      rational& operator-=(rational const& x);
      rational& operator*=(rational const& x);
      rational& operator/=(rational const& x);

      rational& fma(rational const& x, rational const& y);  // Fused multiply-add.  Does *this += x * y;

      int sign() const;

      // this is unfortunate - need some sort of bigint_ref class to avoid the copy.   not yet implemented.
      bigint numerator() const;
      bigint denominator() const;

      std::string to_string() const;

      double to_double() const;

      // IMPLEMENTATION USE ONLY
      mpq_t const& GetRepresentation() const { return Data; }
      mpq_t& GetRepresentation() { return Data; }

   private:
      mpq_t Data;

      // tagged constructors for optimization
      struct abs_tag {};        // absolute value
      struct negate_tag {};     // unary negation
      struct add_tag {};        // binary addition
      struct sub_tag {};        // binary subtraction
      struct mul_tag {};        // binary multiplication
      struct div_tag {};        // binary division

      rational(rational const& x, abs_tag);
      rational(rational const& x, negate_tag);
      rational(rational const& x, rational const& y, add_tag);
      rational(rational const& x, rational const& y, sub_tag);
      rational(rational const& x, rational const& y, mul_tag);
      rational(rational const& x, rational const& y, div_tag);

      rational(bigint const& x, bigint const& y, div_tag);

   // relational operators
   friend bool operator==(rational const& x, rational const& y);
   friend bool operator!=(rational const& x, rational const& y);
   friend bool operator>(rational const& x, rational const& y);
   friend bool operator<(rational const& x, rational const& y);
   friend bool operator>=(rational const& x, rational const& y);
   friend bool operator<=(rational const& x, rational const& y);

   // arithmetic
   friend rational abs(rational const& x);
   friend rational operator-(rational const& x);
   friend rational operator+(rational const& x, rational const& y);
   friend rational operator-(rational const& x, rational const& y);
   friend rational operator*(rational const& x, rational const& y);
   friend rational operator/(rational const& x, rational const& y);

   // integer division -> rational
   friend rational operator/(bigint const& x, bigint const& y);

   friend int sign(rational const& x);

   // IOStreams I/O
#if defined(ENABLE_PERSISTENT_IO)
   friend PStream::opstream& operator<<(PStream::opstream& out, rational const& x);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, rational& x);
#endif

   friend std::ostream& operator<<(std::ostream& out, rational const& x);
   friend std::istream& operator>>(std::istream& in, rational& x);
};

// std::swap

} // namespace

namespace std
{
   template<>
   inline
   void swap<gmp::rational>(gmp::rational& x, gmp::rational& y)
   {
      x.swap(y);
   }
} // namespace std

namespace gmp
{

//
// inlines
//

inline
rational::rational()
{
  mpq_init(Data);
}

inline
rational::rational(rational const& x)
{
   mpq_init(Data);
   mpq_set(Data, x.Data);
}

inline
rational::rational(int x)
{
   mpq_init(Data);
   mpq_set_si(Data, x, 1);
}

inline
rational::rational(long x)
{
   mpq_init(Data);
   mpq_set_si(Data, x, 1);
}

inline
rational::rational(unsigned long x)
{
   mpq_init(Data);
   mpq_set_ui(Data, x, 1);
}

inline
rational::rational(bigint const& x)
{
   mpq_init(Data);
   mpq_set_z(Data, x.GetRepresentation());
}

inline
rational::rational(double x)
{
   mpq_init(Data);
   mpq_set_d(Data, x);
}

// tagged ctors

inline
rational::rational(rational const& x, rational::abs_tag)
{
   mpq_init(Data);
   mpz_abs(mpq_numref(Data), mpq_numref(x.Data));
   mpz_set(mpq_denref(Data), mpq_denref(x.Data));
}

inline
rational::rational(rational const& x, rational::negate_tag)
{
   mpq_init(Data);
   mpz_neg(mpq_numref(Data), mpq_numref(x.Data));
   mpz_set(mpq_denref(Data), mpq_denref(x.Data));
}

inline
rational::rational(rational const& x, rational const& y, rational::add_tag)
{
   mpq_init(Data);
   mpq_add(Data, x.Data, y.Data);
}

inline
rational::rational(rational const& x, rational const& y, rational::sub_tag)
{
   mpq_init(Data);
   mpq_sub(Data, x.Data, y.Data);
}

inline
rational::rational(rational const& x, rational const& y, rational::mul_tag)
{
   mpq_init(Data);
   mpq_mul(Data, x.Data, y.Data);
}

inline
rational::rational(rational const& x, rational const& y, rational::div_tag)
{
   mpq_init(Data);
   mpq_div(Data, x.Data, y.Data);
}

inline
rational::rational(bigint const& x, bigint const& y, rational::div_tag)
{
   mpq_init(Data);
   mpq_set_num(Data, x.GetRepresentation());
   mpq_set_den(Data, y.GetRepresentation());
   mpq_canonicalize(Data);
}

// destructor

inline
rational::~rational()
{
   mpq_clear(Data);
}

// assignment

inline
rational&
rational::operator=(rational const& x)
{
   mpq_set(Data, x.Data);
   return *this;
}

inline
rational&
rational::operator=(long x)
{
   mpq_set_si(Data, x, 1);
   return *this;
}

inline
rational&
rational::operator=(int x)
{
   mpq_set_si(Data, x, 1);
   return *this;
}

inline
rational&
rational::operator=(unsigned long x)
{
   mpq_set_ui(Data, x, 1);
   return *this;
}

// swap

inline
void rational::swap(rational& x)
{
   mpq_swap(Data, x.Data);
}

// extractors

// not yet implemented

// arithmetic

inline
rational&
rational::operator+=(rational const& x)
{
   mpq_add(Data, Data, x.Data);
   return *this;
}

inline
rational&
rational::operator-=(rational const& x)
{
   mpq_sub(Data, Data, x.Data);
   return *this;
}

inline
rational&
rational::operator*=(rational const& x)
{
   mpq_mul(Data, Data, x.Data);
   return *this;
}

inline
rational&
rational::operator/=(rational const& x)
{
   mpq_div(Data, Data, x.Data);
   return *this;
}

inline
int
rational::sign() const
{
   return mpq_sgn(Data);
}

// conversion to string

inline
std::string rational::to_string() const
{
   // get the numerator as a string
   size_t num_size = mpz_sizeinbase(mpq_numref(Data), 10) + 2;
   std::vector<char> num_buffer(num_size);
   mpz_get_str(&num_buffer[0], 10, mpq_numref(Data));

   // get the denominator as a string
   size_t den_size = mpz_sizeinbase(mpq_denref(Data), 10) + 2;
   std::vector<char> den_buffer(den_size);
   mpz_get_str(&den_buffer[0], 10, mpq_denref(Data));

   // put it all together
   return std::string(&num_buffer[0]) + '/' + std::string(&den_buffer[0]);
}

// conversion to double

inline
double rational::to_double() const
{
   return mpq_get_d(Data);
}

//
// non-members
//

// relational operators

inline
bool operator==(rational const& x, rational const& y)
{
   return mpq_equal(x.GetRepresentation(), y.GetRepresentation());
}

inline
bool operator!=(rational const& x, rational const& y)
{
   return mpq_equal(x.GetRepresentation(), y.GetRepresentation()) == 0;
}

inline
bool operator>(rational const& x, rational const& y)
{
   return mpq_cmp(x.GetRepresentation(), y.GetRepresentation()) > 0;
}

inline
bool operator<(rational const& x, rational const& y)
{
   return mpq_cmp(x.GetRepresentation(), y.GetRepresentation()) < 0;
}

inline
bool operator>=(rational const& x, rational const& y)
{
   return mpq_cmp(x.GetRepresentation(), y.GetRepresentation()) >= 0;
}

inline
bool operator<=(rational const& x, rational const& y)
{
   return mpq_cmp(x.GetRepresentation(), y.GetRepresentation()) <= 0;
}

// arithemetic

inline
rational abs(rational const& x)
{
  return rational(x, rational::abs_tag());
}

inline
rational operator-(rational const& x)
{
   return rational(x, rational::negate_tag());
}

inline
rational operator+(rational const& x, rational const& y)
{
   return rational(x,y,rational::add_tag());
}

inline
rational operator-(rational const& x, rational const& y)
{
   return rational(x,y,rational::sub_tag());
}

inline
rational operator*(rational const& x, rational const& y)
{
   return rational(x,y,rational::mul_tag());
}

inline
rational operator/(rational const& x, rational const& y)
{
   return rational(x,y,rational::div_tag());
}

inline
rational operator/(bigint const& x, bigint const& y)
{
   return rational(x,y,rational::div_tag());
}

inline
std::ostream& operator<<(std::ostream& out, rational const& x)
{
   return out << x.to_string();
}

inline
int sign(rational const& x)
{
  return mpq_sgn(x.GetRepresentation());
}

} // namespace gmp

#endif
