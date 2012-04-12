// -*- C++ -*- $Id$

/*
  gmpint.h

  Created 2000-10-14 

  Wraps the GNU-MP library for integers

  It would be nice to eliminate the dependency on gmp.h (lots of global namespace pollution!).
  Can't do that without removing all the inlines though.
*/

#if !defined(GMPINT_H_FHU34T734T785R7Y38734HURHUI54T7YEHIU)
#define GMPINT_H_FHU34T734T785R7Y38734HURHUI54T7YEHIU

#include <algorithm>  // we specialize std::swap
#include <stdexcept>
#include <vector>
#include <stdio.h>
//#if defined(__DECCXX)
using std::FILE;
//#endif
//#include "pstream.h"
#include "gmp.h"
#include "trace.h"

namespace gmp
{

class bigint
{
   public:
      bigint();
      bigint(bigint const& x);
      bigint(int x);
      bigint(long x);
      bigint(unsigned long x);

      bigint(mpz_t const& x);

      explicit bigint(char const* x);
      explicit bigint(std::string const& x);

      ~bigint();

      bigint& operator=(bigint const& x);
      bigint& operator=(long x);
      bigint& operator=(unsigned long x);

      void swap(bigint& x);

      bigint& operator+=(bigint const& x);
      bigint& operator-=(bigint const& x);
      bigint& operator*=(bigint const& x);

      bigint& fma(bigint const& x, bigint const& y);  // Fused multiply-add.  Does *this += x * y;

      // returns true if this number is odd / even
      bool is_odd() const;
      bool is_even() const;

      int sign() const;

      std::string to_string() const;

      // IMPLEMENTATION USE ONLY
      mpz_t const& GetRepresentation() const { return Data; }
      mpz_t& GetRepresentation() { return Data; }

   private:
      mpz_t Data;

      // tagged constructors for optimization
      struct abs_tag {};        // absolute value
      struct negate_tag {};     // unary negation
      struct factorial_tag {};  // factorial
      struct binomial_tag {};   // binomial (n k)
      struct add_tag {};        // binary addition
      struct sub_tag {};        // binary subtraction
      struct mul_tag {};        // binary multiplication

      bigint(bigint const& x, abs_tag);
      bigint(bigint const& x, negate_tag);
      bigint(unsigned long x, factorial_tag);
      bigint(unsigned long n, unsigned long k, binomial_tag);
      bigint(bigint const& n, unsigned long k, binomial_tag);
      bigint(bigint const& x, bigint const& y, add_tag);
      bigint(bigint const& x, bigint const& y, sub_tag);
      bigint(bigint const& x, bigint const& y, mul_tag);

   // relational operators
   friend bool operator==(bigint const& x, bigint const& y);
   friend bool operator>(bigint const& x, bigint const& y);
   friend bool operator<(bigint const& x, bigint const& y);
   friend bool operator>=(bigint const& x, bigint const& y);
   friend bool operator<=(bigint const& x, bigint const& y);

   friend bigint abs(bigint const& x);
   friend bigint operator-(bigint const& x);

   // arithmetic binary operators
   friend bigint operator+(bigint const& x, bigint const& y);
   friend bigint operator-(bigint const& x, bigint const& y);
   friend bigint operator*(bigint const& x, bigint const& y);

   // division not yet implemented - what notation to use?

   // factorial operations.  all integer types are included here to avoid possible
   // overloading problems.
   friend bigint factorial(unsigned long x);
   friend bigint factorial(long x);
   friend bigint factorial(unsigned int x);
   friend bigint factorial(int x);
   friend bigint factorial(unsigned short x);
   friend bigint factorial(short x);

   friend bigint binomial(unsigned long n, unsigned long k);
   friend bigint binomial(bigint const& n, unsigned long k);

   // persistent I/O
#if defined(DEFINE_PERSISTENT_IO)
   friend PStream::opstream& operator<<(PStream::opstream& out, bigint const& x);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, bigint& x);
#endif

   // IOStreams I/O
   friend std::ostream& operator<<(std::ostream& out, bigint const& x);
   friend std::istream& operator>>(std::istream& in, bigint& x);

   // NOTE: the whitespace delimiting is different for the stream extractor and the string constructor.
   // the string ctor uses GMP semantics, but the extractor uses standard C++ string extractor semantics.
};

} // namespace gmp

namespace std
{
   template<>
   inline
   void swap<gmp::bigint>(gmp::bigint& x, gmp::bigint& y)
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
bigint::bigint()
{
  mpz_init(Data);
}

inline
bigint::bigint(bigint const& x)
{
   mpz_init_set(Data, x.Data);
}

inline
bigint::bigint(int x)
{
   mpz_init_set_si(Data, x);
}

inline
bigint::bigint(long x)
{
   mpz_init_set_si(Data, x);
}

inline
bigint::bigint(unsigned long x)
{
   mpz_init_set_ui(Data, x);
}

inline
bigint::bigint(mpz_t const& x)
{
   mpz_init_set(Data, x);
}

inline
bigint::bigint(char const* x)
{
   mpz_init_set_str(Data, x, 10);
}

inline
bigint::bigint(std::string const& x)
{
   mpz_init_set_str(Data, x.c_str(), 10);
}

// tagged ctors

inline
bigint::bigint(bigint const& x, bigint::abs_tag)
{
   mpz_init(Data);
   mpz_abs(Data, x.Data);
}

inline
bigint::bigint(bigint const& x, bigint::negate_tag)
{
   mpz_init(Data);
   mpz_neg(Data, x.Data);
}

inline
bigint::bigint(unsigned long x, bigint::factorial_tag)
{
   mpz_init(Data);
   mpz_fac_ui(Data, x);
}

inline
bigint::bigint(unsigned long n, unsigned long k, bigint::binomial_tag)
{
   mpz_init(Data);
   mpz_bin_uiui(Data, n, k);
}

inline
bigint::bigint(bigint const& n, unsigned long k, bigint::binomial_tag)
{
   mpz_init(Data);
   mpz_bin_ui(Data, n.Data, k);
}

inline
bigint::bigint(bigint const& x, bigint const& y, bigint::add_tag)
{
   mpz_init(Data);
   mpz_add(Data, x.Data, y.Data);
}

inline
bigint::bigint(bigint const& x, bigint const& y, bigint::sub_tag)
{
   mpz_init(Data);
   mpz_sub(Data, x.Data, y.Data);
}

inline
bigint::bigint(bigint const& x, bigint const& y, bigint::mul_tag)
{
   mpz_init(Data);
   mpz_mul(Data, x.Data, y.Data);
}

// destructor

inline
bigint::~bigint()
{
   mpz_clear(Data);
}

// assignment

inline
bigint&
bigint::operator=(bigint const& x)
{
   mpz_set(Data, x.Data);
   return *this;
}

inline
bigint&
bigint::operator=(long x)
{
   mpz_set_si(Data, x);
   return *this;
}

inline
bigint&
bigint::operator=(unsigned long x)
{
   mpz_set_ui(Data, x);
   return *this;
}

// swap

inline
void bigint::swap(bigint& x)
{
   mpz_swap(Data, x.Data);
}

// arithmetic

inline
bigint& 
bigint::operator+=(bigint const& x)
{
   mpz_add(Data, Data, x.Data);  // we are allowed to alias parameters in GNU MP
   return *this;
}

inline
bigint& 
bigint::operator-=(bigint const& x)
{
   mpz_sub(Data, Data, x.Data); 
   return *this;
}

inline
bigint& 
bigint::operator*=(bigint const& x)
{
   mpz_mul(Data, Data, x.Data);  
   return *this;
}

inline
bool 
bigint::is_odd() const
{
   return mpz_odd_p(Data);
}

inline
bool 
bigint::is_even() const
{
   return mpz_even_p(Data);
}

inline
int
bigint::sign() const
{
   return mpz_sgn(Data);
}

// conversion to string

inline
std::string bigint::to_string() const
{
   size_t size = mpz_sizeinbase(Data, 10) + 2;
   std::vector<char> buffer(size);
   mpz_get_str(&buffer[0], 10, Data);
   return std::string(&buffer[0]);
}

// non-members

// comparision

inline
bool operator==(bigint const& x, bigint const& y)
{
   return mpz_cmp(x.Data, y.Data) == 0;
}

inline
bool operator>(bigint const& x, bigint const& y)
{
   return mpz_cmp(x.Data, y.Data) > 0;
}

inline
bool operator<(bigint const& x, bigint const& y)
{
   return mpz_cmp(x.Data, y.Data) < 0;
}

inline
bool operator>=(bigint const& x, bigint const& y)
{
   return mpz_cmp(x.Data, y.Data) >= 0;
}

inline
bool operator<=(bigint const& x, bigint const& y)
{
   return mpz_cmp(x.Data, y.Data) <= 0;
}

// arithmetic

inline
bigint
abs(bigint const& x)
{
  return bigint(x, bigint::abs_tag());
}

inline
bigint
operator-(bigint const& x)
{
   return bigint(x, bigint::negate_tag());
}

inline
bigint operator+(bigint const& x, bigint const& y)
{
   return bigint(x,y,bigint::add_tag());
}

inline
bigint operator-(bigint const& x, bigint const& y)
{
   return bigint(x,y,bigint::sub_tag());
}

inline
bigint operator*(bigint const& x, bigint const& y)
{
   return bigint(x,y,bigint::mul_tag());
}

// factorial

inline
bigint 
factorial(unsigned long x)
{
   return bigint(x, bigint::factorial_tag());
}

inline
bigint 
factorial(long x)
{
   DEBUG_PRECONDITION(x >= 0);
   return bigint(x, bigint::factorial_tag());
}

inline
bigint 
factorial(unsigned int x)
{
   return bigint(x, bigint::factorial_tag());
}

inline
bigint 
factorial(int x)
{
   DEBUG_PRECONDITION(x >= 0);
   return bigint(x, bigint::factorial_tag());
}

inline
bigint 
factorial(unsigned short x)
{
   return bigint(x, bigint::factorial_tag());
}

inline
bigint 
factorial(short x)
{
   DEBUG_PRECONDITION(x >= 0);
   return bigint(x, bigint::factorial_tag());
}

// binomial

inline
bigint
binomial(unsigned long n, unsigned long k)
{
   return bigint(n, k, bigint::binomial_tag());
}

inline
bigint
binomial(bigint const& n, unsigned long k)
{
   return bigint(n, k, bigint::binomial_tag());
}

// stream I/O

inline
std::ostream& operator<<(std::ostream& out, bigint const& x)
{
   return out << x.to_string();
}

inline
std::istream& operator>>(std::istream& in, bigint& x)
{
   std::string s;
   in >> s;
   x = bigint(s);
   return in;
}

} // namespace gmp

#endif
