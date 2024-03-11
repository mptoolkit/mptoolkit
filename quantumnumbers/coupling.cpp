// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/coupling.cpp
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
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

#include "coupling.h"
#include "common/gmprational.h"
#include "common/sortsearch.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <unordered_map>
#include <mutex>
using gmp::rational;
using gmp::bigint;
using gmp::factorial;

/* See: [Thompson, Atlas for Computing Mathematical Functions] */

rational
delta(half_int ta, half_int tb, half_int tc)
{
   bigint f1 = gmp::factorial(to_int(ta+tb-tc));
   bigint f2 = gmp::factorial(to_int(ta+tc-tb));
   bigint f3 = gmp::factorial(to_int(tb+tc-ta));
   bigint f4 = gmp::factorial(to_int(ta+tb+tc)+1);

   return f1 * f2 * f3 / f4;
}

bool check_m(half_int ja, half_int jb, half_int jc, half_int ma, half_int mb, half_int mc)
{
   return (abs(ma) <= ja) && (abs(mb) <= jb) && (abs(mc) <= jc)
     && is_integral(ja+ma) && is_integral(jb+mb) && is_integral(jc+mc);
}

// low level Racah coefficient
double Racah_NoCache(half_int ja, half_int jb, half_int jc, half_int jd, half_int je, half_int jf)
{
   double norm = sqrt((delta(ja,jb,je)*delta(jc,jd,je)*delta(ja,jc,jf)*delta(jb,jd,jf)).to_double());

   half_int zmin = max4(ja+jb+je, jc+jd+je, ja+jc+jf, jb+jd+jf);
   half_int zmax = min3(ja+jb+jc+jd, ja+jd+je+jf, jb+jc+je+jf);

   rational sum_pos = 0, sum_neg = 0;

   // this calculation splits up the +ve and -ve parts.  There is no point doing this
   // now that we are using rational arithmetic.
   int phase = minus1pow(to_int(ja+jb+jc+jd+zmin));
   for (half_int z = zmin; z <= zmax; ++z)
   {
      bigint n = factorial(to_int(z+1));

      bigint d1 = factorial(to_int(z-ja-jb-je));
      bigint d2 = factorial(to_int(z-jc-jd-je));
      bigint d3 = factorial(to_int(z-ja-jc-jf));
      bigint d4 = factorial(to_int(z-jb-jd-jf));

      bigint d5 = factorial(to_int(ja+jb+jc+jd-z));
      bigint d6 = factorial(to_int(ja+jd+je+jf-z));
      bigint d7 = factorial(to_int(jb+jc+je+jf-z));

      rational term = n / (d1*d2*d3*d4*d5*d6*d7);

      if (phase > 0) sum_pos += term;
      else sum_neg += term;

      phase = -phase;
   }

   return norm * (sum_pos - sum_neg).to_double();
}

// the actual Racah function is implemented in terms of Coupling6j
double Racah(half_int ja, half_int jb, half_int jc, half_int jd, half_int je, half_int jf)
{
   half_int phase = ja + jb + jc + jd;
   if (!phase.is_integral()) return 0;
   return minus1pow(to_int(phase)) * Coupling6j(ja,jb,je,jd,jc,jf);
}

// struct to hold the 6j labels in a hash table
struct Coefficients6j
{
   int j1, j2, j3, j4, j5, j6;

   Coefficients6j(half_int j1_, half_int j2_, half_int j3_, half_int j4_, half_int j5_, half_int j6_);

   // transform to a canonical form, making use of the symmetries.
   // If the coefficients were already in canonical form, we return true.
   bool Canonicalize();
};

inline
bool operator==(Coefficients6j const& x, Coefficients6j const& y)
{
   return x.j1 == y.j1 && x.j2 == y.j2 && x.j3 == y.j3 &&
     x.j4 == y.j4 && x.j5 == y.j5 && x.j6 == y.j6;
}

inline
Coefficients6j::Coefficients6j(half_int j1_, half_int j2_, half_int j3_,
                               half_int j4_, half_int j5_, half_int j6_)
   : j1(j1_.twice()), j2(j2_.twice()), j3(j3_.twice()), j4(j4_.twice()), j5(j5_.twice()), j6(j6_.twice())
{
}

bool Coefficients6j::Canonicalize()
{
   // Use the symmetries of the 6j coefficients to get a 'cannonical'
   // form.  the symmetry group of the 6j's is of order 24.  (Other symmetries mix/change the coefficients.)
   // We can swap pairs (j1,j4), (j2,j5), (j3,j6) arbitarily, and we can swap any two sets of
   // (j1 <--> j4), (j2 <--> j5), (j3 <--> j6).

   using std::swap; using std::min; using std::max;

   bool ColumnMix = true, RowMix = true;

   // firstly, we can freely swap columns (j1,j4), (j2,j5) and (j3,j6) among themselves.
   // Put them in decreasing order by maximum
   int max14 = max(j1,j4);
   int max25 = max(j2,j5);
   int max36 = max(j3,j6);

   int min14 = min(j1,j4);
   int min25 = min(j2,j5);
   int min36 = min(j3,j6);

   // put them in order 14 >= 25 >= 36.  We don't care about the values of the min/max after this
   // sequence, so don't bother swapping them
   if (max25 > max14 || (max25 == max14 && min25 > min14))
   {
      // 25 > 14
      if (max36 > max25 || (max36 == max25 && min36 > min25))
      {
         // 36 > 25 > 14.  swap(36,14)
         swap(j1, j3);
         swap(j4, j6);
         // swap(max14, max36); - not needed
         // swap(min14, min36); - not needed
      }
      else if (max36 > max14 || (max36 == max14 && min36 > min14))
      {
         // 25 >= 36 > 14.  rotate.
         int w = j1, x = j4;
         j1 = j2;
         j4 = j5;

         j2 = j3;
         j5 = j6;

         j3 = w;
         j6 = x;
      }
      else
      {
         // 25 > 14 >= 36.  swap(25,14)
         swap(j1, j2);
         swap(j4, j5);
      }
   }
   else
   {
      // 14 >= 25
      if (max36 > max14 || (max36 == max14 && min36 > min14))
      {
         // 36 > 14 >= 25.  rotate
         int w = j2, x = j5;

         j2 = j1;
         j5 = j4;

         j1 = j3;
         j4 = j6;

         j3 = w;
         j6 = x;
      }
      else if (max36 > max25 || (max36 == max25 && min36 > min25))
      {
         // 14 >= 36 > 25.  swap(36, 25)
         swap(j3, j2);
         swap(j6, j5);
      }
      else
      {
         // final combination is 14 >= 25 >= 36.  Already ordered, no action.
         ColumnMix = false;
      }
   }

   // Now we can swap within pairs.  If we do make a swap, we must swap exactly two pairs.
   // The strategy here is to put the largest elements on top.  If we swap (j1,j4) then
   // we also need to swap one of (j2,j5) or (j3,j6).  If we do not swap (j1,j4) then
   // we decide whether to swap (j2,j5) (and hence also (j3,j6)).
   // An interesting case is when one of the columns have equal entries;
   // here we can freely swap entries in any column.

   if (j1 < j4)
   {
      swap(j1, j4);

      // now decide which of (j2,j5) or (j3,j6) to swap
      if (j2 < j5)
      {
         swap(j2, j5);
      }
      else if (j3 < j6 || (j5 < j2 && j6 < j3))
      {
         swap(j3, j6);
      }
      else
      {
         // in this case, one column has equal entries so we don't need to do anything.
         CHECK(j2 == j5 || j3 == j6);
      }
   }
   else if (j2 < j5)
   {
      swap(j2,j5);
      if (j1 != j4)
         swap(j3,j6);
   }
   else if (j3 < j6 && (j1 == j4 || j2 == j5))
   {
      swap(j3, j6);
   }
   else
   {
      // if we get here then we don't need to swap any pairs
      RowMix = false;
   }

   return !ColumnMix && !RowMix;
}

double
Coupling6j_NoCache(Coefficients6j const& x)
{
   int sign = minus1pow((x.j1 + x.j2 + x.j4 + x.j5) / 2);
   return sign * Racah_NoCache(from_twice(x.j1), from_twice(x.j2),
                               from_twice(x.j5), from_twice(x.j4),
                               from_twice(x.j3), from_twice(x.j6));
}

std::ostream& operator<<(std::ostream& out, Coefficients6j const& x) // for debugging
{
   return out << '{' << x.j1 << ',' << x.j2 << ',' << x.j3 << ','
              << x.j4 << ',' << x.j5 << ',' << x.j6 << '}';
}

// specialization of the ext::hash<> function for Coefficients6j

namespace std
{

template <>
struct hash<Coefficients6j>
{
   static int const HashBits6j = std::numeric_limits<size_t>::digits / 6;

   size_t operator()(Coefficients6j const& c) const
   {
      size_t ret = c.j1;
      ret <<= HashBits6j;
      ret += c.j2;
      ret <<= HashBits6j;
      ret += c.j3;
      ret <<= HashBits6j;
      ret += c.j4;
      ret <<= HashBits6j;
      ret += c.j5;
      ret <<= HashBits6j;
      ret += c.j6;
      return ret;
   }
};

} // namespace std

typedef std::unordered_map<Coefficients6j, double> Hash6jType;
Hash6jType HashTable;
std::mutex HashMutex;

double Coupling6j(half_int j1, half_int j2, half_int j3, half_int j4, half_int j5, half_int j6)
{

   // do some checking of parameters
   if (j1 < 0 || j2 < 0 || j3 < 0 || j4 < 0 || j5 < 0 || j6 < 0) throw std::domain_error("negative j in Coupling6j()");

   if (!is_triangle(j1,j2,j3) || !is_triangle(j4,j5,j3) || !is_triangle(j1,j5,j6) || !is_triangle(j2,j4,j6))
     return 0;

   // hashing is probably faster than canonicalizing the 6j's, so we look up the hash table first
   // for the non-canonical form.
   Coefficients6j Coeff(j1, j2, j3, j4, j5, j6);
#if defined(CONFIG_CANONICALIZE_6J)
   Coeff.Canonicalize();
#endif

   std::lock_guard<std::mutex> MyLock(HashMutex);

   Hash6jType::iterator I = HashTable.find(Coeff);
   if (I == HashTable.end())
   {
#if !defined(CONFIG_CANONICALIZE_6J)
      // not in the hash table.  Canonicalize the coefficients.
      Coefficients6j CanonicalCoeff(Coeff);
      bool WasCanonical = CanonicalCoeff.Canonicalize();
      if (!WasCanonical)
      {
         I = HashTable.find(CanonicalCoeff);
         if (I == HashTable.end())
         {
            // neither the original nor the canonical entry was in the table
            double x = Coupling6j_NoCache(CanonicalCoeff);
            HashTable[CanonicalCoeff] = x;
            HashTable[Coeff] = x;
            return x;
         }
         else
         {
            // the canonical form is in the table
            double x = I->second;
            HashTable[Coeff] = x;
            return x;
         }
      }
      else
#endif
      {
         // it was already in canonical form
         double x = Coupling6j_NoCache(Coeff);
         HashTable[Coeff] = x;
         return x;
      }
   }
   // else the coefficient was already in the table, no canonicalization needed
   return I->second;
}

double Coupling9j(half_int j11, half_int j12, half_int j13,
                  half_int j21, half_int j22, half_int j23,
                  half_int j31, half_int j32, half_int j33)
{
  //   std::cout<<std::setw(3)<<j11<<' '<<std::setw(3)<<j12<<' '<<std::setw(3)<<j13<<'\n'
  //            <<std::setw(3)<<j21<<' '<<std::setw(3)<<j22<<' '<<std::setw(3)<<j23<<'\n'
  //            <<std::setw(3)<<j31<<' '<<std::setw(3)<<j32<<' '<<std::setw(3)<<j33<<std::endl;

   if (j11 < 0 || j12 < 0 || j13 < 0 ||
       j21 < 0 || j22 < 0 || j23 < 0 ||
       j31 < 0 || j32 < 0 || j33 < 0)
   {
      throw std::domain_error("negative j in Coupling9j()");
   }

   if (!is_triangle(j11, j12, j13) ||
       !is_triangle(j21, j22, j23) ||
       !is_triangle(j31, j32, j33) ||
       !is_triangle(j11, j21, j31) ||
       !is_triangle(j12, j22, j32) ||
       !is_triangle(j13, j23, j33))
   {
     //      std::cout << 0 << std::endl;
      return 0;
   }

   half_int kmin = max3(abs(j21-j32), abs(j12-j23), abs(j11-j33));
   half_int kmax = min3(j21+j32, j12+j23, j11+j33);

   double sum_neg = 0, sum_pos = 0;
   int phase = minus1pow(kmin.twice());

   for (half_int k = kmin; k <= kmax; k += 1)
   {
      double s1 = Coupling6j(j11,j21,j31, j32,j33,k);
      double s2 = Coupling6j(j12,j22,j32, j21,k,j23);
      double s3 = Coupling6j(j13,j23,j33, k,j11,j12);

      double term = s1 * s2 * s3;

      if(term > 0)
      {
        sum_pos += (k.twice() + 1) * term;
      }
      else
      {
        sum_neg += (k.twice() + 1) * term;
      }
   }

   double Result = sum_pos + sum_neg;
   //TRACE_IF(std::abs(Result / (sum_pos-sum_neg)) < std::numeric_limits<double>::epsilon()*10)("numerical zero 9j")(Result);
   if (std::abs(Result / (sum_pos-sum_neg)) < std::numeric_limits<double>::epsilon()*10)
      return 0.0;
   // else
   return phase * (sum_pos + sum_neg);
}

namespace
{

rational CGSummation(half_int j1, half_int m1, half_int j2, half_int m2, half_int j)
{
   long KMax = to_int(std::min(std::min(j1+j2-j, j1-m1), j2+m2));
   long KMin = to_int(std::max(half_int(0), std::max(j2-j-m1, j1-j+m2)));

   rational Sum = 0;

   long j1j2_j = to_int(j1 + j2 - j);
   long j1_m1 = to_int(j1 - m1);
   long j2m2 = to_int(j2 + m2);
   long j_j2m1 = to_int(j - j2 + m1);
   long j_j1_m2 = to_int(j - j1 - m2);

   for (long k = KMin; k <= KMax; ++k)
   {
      Sum += minus1pow(k) / (gmp::factorial(k) * gmp::factorial(j1j2_j - k) * gmp::factorial(j1_m1 - k)
        * gmp::factorial(j2m2 - k) * gmp::factorial(j_j2m1 + k) * gmp::factorial(j_j1_m2 + k));
   };

   return Sum;
}

double PreFactor(
half_int j1, half_int m1,
half_int j2, half_int m2,
half_int j, half_int m)
{
   rational T = to_int(2*j+1) * gmp::factorial(to_int(j1+j2-j)) * gmp::factorial(to_int(j1-j2+j)) *
     gmp::factorial(to_int(-j1+j2+j))
               / gmp::factorial(to_int(j1+j2+j+1));

   return sqrt((T * gmp::factorial(to_int(j1+m1))*gmp::factorial(to_int(j1-m1))*gmp::factorial(to_int(j2+m2))*
                          gmp::factorial(to_int(j2-m2))*
                          gmp::factorial(to_int(j+m))*gmp::factorial(to_int(j-m))).to_double());
}

// returns the square of the prefactor as a rational
rational PreFactorSq(
half_int j1, half_int m1,
half_int j2, half_int m2,
half_int j, half_int m)
{
   rational T = to_int(2*j+1) * gmp::factorial(to_int(j1+j2-j)) * gmp::factorial(to_int(j1-j2+j)) *
     gmp::factorial(to_int(-j1+j2+j))
               / gmp::factorial(to_int(j1+j2+j+1));

   return T * gmp::factorial(to_int(j1+m1))*gmp::factorial(to_int(j1-m1))*gmp::factorial(to_int(j2+m2))*
                          gmp::factorial(to_int(j2-m2))*
      gmp::factorial(to_int(j+m))*gmp::factorial(to_int(j-m));
}

} // namespace

double ClebschGordan(
half_int j1, half_int m1,
half_int j2, half_int m2,
half_int j,  half_int m)
{
   // first, the paranoid check
   if (j1 < 0 || j2 < 0 || j < 0)
      throw std::invalid_argument("negative angular momenta encountered in function ClebschGordan");

   // any coefficients outside the allowable ranges are zero
   if (m1 + m2 != m || j < abs(j1 - j2) || j > j1 + j2 || abs(m1) > j1 || abs(m2) > j2
                    || abs(m) > j) return 0;

   // Also, the half-integer parts must match
   if (!is_integral(j1 + m1) || !is_integral(j2+m2) || !is_integral(j+m)
        || !is_integral(j1+j2-j)) return 0;

   if (j == 0)
   {
      // assert(j1 == j2 && m1 == -m2);
      return minus1pow(to_int(j1-m1)) / sqrt(double(j1.twice() + 1));
   };
   // else

   return PreFactor(j1,m1,j2,m2,j,m) * CGSummation(j1,m1,j2,m2,j).to_double();
}

std::pair<rational, rational>
ClebschGordanSquared(half_int j1, half_int m1,
                     half_int j2, half_int m2,
                     half_int j,  half_int m)
{
   // first, the paranoid check
   if (j1 < 0 || j2 < 0 || j < 0)
      throw std::invalid_argument("negative angular momenta encountered in function ClebschGordanSquared");

   // any coefficients outside the allowable ranges are zero
   if (m1 + m2 != m || j < abs(j1 - j2) || j > j1 + j2 || abs(m1) > j1 || abs(m2) > j2
                    || abs(m) > j) return std::make_pair(0,0);

   // Also, the half-integer parts must match
   if (!is_integral(j1 + m1) || !is_integral(j2+m2) || !is_integral(j+m)
        || !is_integral(j1+j2-j)) return std::make_pair(0,0);

   return std::make_pair(CGSummation(j1,m1,j2,m2,j), PreFactorSq(j1,m1,j2,m2,j,m));
}
