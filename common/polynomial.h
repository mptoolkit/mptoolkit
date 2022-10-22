// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/polynomial.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// A simple polynomial class over a templated coefficient field.
// The domain is not specified in this class, rather, this class defines
// an abstract polynomial of the form
// sum_i a_i X^i, where X is some abstract symbol.
// There are templated functions to evaluate the polynomial and its derivatives.
//
// Addition and subtraction are defined, as addition of the coefficients
// of each degree.

#if !defined(MPTOOLKIT_COMMON_POLYNOMIAL_H)
#define MPTOOLKIT_COMMON_POLYNOMIAL_H

#include <map>
#include <ostream>
#include <cmath>

template <typename CoefficientField>
class Polynomial
{
   public:
      typedef CoefficientField coefficient_type;
      typedef coefficient_type data_type;
      typedef int key_type;

      typedef std::map<int, coefficient_type> container_type;

      typedef typename container_type::iterator       iterator;
      typedef typename container_type::const_iterator const_iterator;
      typedef typename container_type::value_type value_type;

      // default constructor makes a zero polynomial
      Polynomial() {}

      // initialize as a constant polynomial
      explicit Polynomial(coefficient_type const& c);

      // returns the coefficient of the term x^n, or a default constructed coefficient
      // if the term is zero
      coefficient_type coefficient(int n) const;

      // returns true if this Polynomial has a term at degree n
      bool has_term(int n) const { return data_.find(n) != data_.end(); }

      // bracket operator P[n] gives the coefficient of the term x^n,
      // as an lvalue.
      coefficient_type& operator[](int n) { return data_[n]; }

      coefficient_type operator[](int n) const { return this->coefficient(n); }

      // evaluate the polynomial at point x.  Returns sum_n coefficient(n) * pow(x,n)
      template <typename T>
      auto evaluate(T const& x) const -> decltype(coefficient_type() * x);

      // equivalent to evaluate()
      template <typename T>
      auto operator()(T const& x) const -> decltype(coefficient_type() * x);

      // evalutes the first derivative
      template <typename T>
      auto derivative(T const& x) const -> decltype(coefficient_type() * x);

      // evaluates the n'th derivative.  n must be non-negative.
      template <typename T>
      auto derivative(T const& x, int n) const -> decltype(coefficient_type() * x);

      bool empty() const { return data_.empty(); }
      int degree() const;

      // returns the number of non-zero coefficients
      int size() const { return int(data_.size()); }

      iterator begin() { return data_.begin(); }
      iterator end() { return data_.end(); }

      const_iterator begin() const { return data_.begin(); }
      const_iterator end() const { return data_.end(); }

      std::map<int, coefficient_type> data_;
};

template <typename CF>
Polynomial<CF>::Polynomial(coefficient_type const& c)
{
   data_[0] = c;
}

template <typename CF>
int Polynomial<CF>::degree() const
{
   if (data_.empty())
      return 0;
   // else
   const_iterator i = data_.end();
   --i;
   return i->first;
}

template <typename CF>
CF Polynomial<CF>::coefficient(int n) const
{
   const_iterator I = data_.find(n);
   if (I == data_.end())
      return coefficient_type();
   return I->second;
}

template <typename CF>
template <typename T>
auto
Polynomial<CF>::evaluate(T const& x) const -> decltype(CF() * x)
{
   using std::pow;
   const_iterator i = this->begin();
   if (i == this->end())
   {
      return decltype(coefficient_type() * x){};
   }
   auto v = i->second * pow(x, i->first);
   ++i;
   while (i != this->end())
   {
      v += i->second * pow(x, i->first);
      ++i;
   }
   return v;
}

template <typename CF>
template <typename T>
auto
Polynomial<CF>::operator()(T const& x) const -> decltype(CF() * x)
{
   return this->evaluate(x);
}

// returns the first r terms of n!, i.e. n! / (n-r)!
inline
double factorial_r(int n, int r)
{
   double v = 1;
   for (int i = n; i > n-r; --i)
   {
      v *= i;
   }
   return v;
}

template <typename CF>
template <typename T>
auto
Polynomial<CF>::derivative(T const& x, int n) const -> decltype(CF() * x)
{
   using std::pow;
   const_iterator i = this->begin();
   // find the first non-zero component
   while (i != this->end() && i->first < n)
      ++i;
   if (i == this->end())
   {
      return decltype(coefficient_type() * x){};
   }
   auto v = i->second * pow(x, i->first - n) * factorial_r(i->first, n);
   ++i;
   while (i != this->end())
   {
      v += i->second * pow(x, i->first - n) * factorial_r(i->first, n);
      ++i;
   }
   return v;
}

template <typename CF>
template <typename T>
auto
Polynomial<CF>::derivative(T const& x) const -> decltype(CF() * x)
{
   return this->derivative(x, 1);
}

template <typename CF>
std::ostream&
operator<<(std::ostream& out, Polynomial<CF> const& x)
{
   if (x.empty())
   {
      return out << "zero polynomial\n";
   }

   for (typename Polynomial<CF>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      out << "Coefficient of x^" << I->first << " = " << I->second << '\n';
   }
   return out;
}

template <typename CF>
Polynomial<CF>&
operator+=(Polynomial<CF>& Poly, Polynomial<CF> const& x)
{
   for (typename Polynomial<CF>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      Poly[I->first] += I->second;
   }
   return Poly;
}

template <typename CF>
Polynomial<CF>&
operator-=(Polynomial<CF>& Poly, Polynomial<CF> const& x)
{
   for (typename Polynomial<CF>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      Poly[I->first] -= I->second;
   }
   return Poly;
}

template <typename CF>
Polynomial<CF>
operator+(Polynomial<CF> const& x, Polynomial<CF> const& y)
{
   Polynomial<CF> Result(x);
   Result += y;
   return Result;
}

template <typename CF>
Polynomial<CF>
operator-(Polynomial<CF> const& x, Polynomial<CF> const& y)
{
   Polynomial<CF> Result(x);
   Result -= y;
   return Result;
}

template <typename CF, typename Scalar>
Polynomial<CF>&
operator*=(Polynomial<CF>& Poly, Scalar x)
{
   for (auto& m : Poly)
   {
      m.second *= x;
   }
   return Poly;
}

template <typename CF, typename Scalar>
Polynomial<CF>
operator*(Polynomial<CF> Poly, Scalar x)
{
   for (auto& m : Poly)
   {
      m.second = m.second * x;
   }
   return Poly;
}

template <typename CF, typename Scalar>
Polynomial<CF>
operator*(Scalar x, Polynomial<CF> Poly)
{
   for (auto& m : Poly)
   {
      m.second = x * m.second;
   }
   return Poly;
}

#endif
