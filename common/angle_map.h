// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/angle_map.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// angle_map.h
// a std::map for angular quantities

#if !defined(MPTOOLKIT_COMMON_ANGLE_MAP_H)
#define MPTOOLKIT_COMMON_ANGLE_MAP_H

#include <complex>
#include <map>
#include <cmath>
#include <ostream>
#include <vector>
#include "common/math_const.h"
#include "common/trace.h"

double const default_angle_resolution = 1E-10;

// Comparitor for complex numbers.  This is so that we can put them in a map,
// the choice of comparison operation is arbitrary
struct CompareComplex
{
   typedef std::complex<double> first_argument_type;
   typedef std::complex<double> second_argument_type;
   typedef bool result_type;
   bool operator()(std::complex<double> const& x, std::complex<double> const& y) const
   {
      return (x.real() < y.real()) || (x.real() == y.real() && x.imag() < y.imag());
   }
};

template <typename T>
class angle_map
{
   private:
      typedef std::map<std::complex<double>, T, CompareComplex> map_type;

   public:
      typedef typename map_type::iterator       iterator;
      typedef typename map_type::const_iterator const_iterator;

      angle_map();
      angle_map(angle_map const& other);

      T& operator[](std::complex<double> const& c);

      T lookup_or_default(std::complex<double> const& c);

      bool has_element(std::complex<double> const& c);

      iterator begin() { return Map.begin(); }
      iterator end() { return Map.end(); }

      const_iterator begin() const { return Map.begin(); }
      const_iterator end() const { return Map.end(); }

   private:
      map_type Map;
      double Resolution;
};

template <typename T>
angle_map<T>::angle_map()
   : Resolution(default_angle_resolution)
{
}


template <typename T>
angle_map<T>::angle_map(angle_map<T> const& Other)
   : Map(Other.Map), Resolution(default_angle_resolution)
{
}

template <typename T>
T&
angle_map<T>::operator[](std::complex<double> const& c)
{
   double r = c.real()*c.real() + c.imag()*c.imag();
   CHECK(std::abs(r - 1.0) <= Resolution)("complex angle is not sufficiently close to the unit circle!")(r);

   // Firstly, we treat angle 0 and pi as special cases
   double angle = atan2(c.imag(), c.real());
   if (-Resolution < angle && angle < Resolution)
   {
      // we have a zero angle
      return Map[std::complex<double>(1.0, 0.0)];
   }
   // else

   if (math_const::pi-Resolution < angle || angle < -math_const::pi+Resolution)
   {
      // we have a pi angle
      return Map[std::complex<double>(-1.0, 0.0)];
   }
   // else

   // Now search through the map
   iterator I = Map.begin();
   while (I != Map.end() && std::abs(c - I->first) > Resolution*Resolution)
   {
      ++I;
   }

   if (I == Map.end())
   {
      // new angle
      DEBUG_TRACE("Adding new angle")(c);
      return Map[c];
   }
   // else
   return I->second;
}

template <typename T>
T
angle_map<T>::lookup_or_default(std::complex<double> const& c)
{
   double r = c.real()*c.real() + c.imag()*c.imag();
   CHECK(std::abs(r - 1.0) <= Resolution)("complex angle is not sufficiently close to the unit circle!")(r);

   // Firstly, we treat angle 0 and pi as special cases
   double angle = atan2(c.imag(), c.real());
   if (-Resolution < angle && angle < Resolution)
   {
      // we have a zero angle
      const_iterator I = Map.find(std::complex<double>(1.0, 0.0));
      return I == Map.end() ? T() : I->second;
   }
   // else

   if (math_const::pi-Resolution < angle || angle < -math_const::pi+Resolution)
   {
      // we have a pi angle
      const_iterator I = Map.find(std::complex<double>(-1.0, 0.0));
      return I == Map.end() ? T() : I->second;
   }
   // else

   // Now search through the map
   iterator I = Map.begin();
   while (I != Map.end() && std::abs(c - I->first) > Resolution*Resolution)
   {
      ++I;
   }

   if (I == Map.end())
   {
      // not found
      return T();
   }
   // else
   return I->second;
}

template <typename T>
bool
angle_map<T>::has_element(std::complex<double> const& c)
{
   double r = c.real()*c.real() + c.imag()*c.imag();
   CHECK(std::abs(r - 1.0) <= Resolution)("complex angle is not sufficiently close to the unit circle!")(r);

   // Firstly, we treat angle 0 and pi as special cases
   double angle = atan2(c.imag(), c.real());
   if (-Resolution < angle && angle < Resolution)
   {
      // we have a zero angle
      return Map.find(std::complex<double>(1.0,0.0)) != Map.end();
   }
   // else

   if (math_const::pi-Resolution < angle || angle < -math_const::pi+Resolution)
   {
      // we have a pi angle
      return Map.find(std::complex<double>(-1.0,0.0)) != Map.end();
   }
   // else

   // Now search through the map
   iterator I = Map.begin();
   while (I != Map.end() && std::abs(c - I->first) > Resolution*Resolution)
   {
      ++I;
   }

   if (I == Map.end())
   {
      return false;
   }
   return true;
}

template <typename T>
inline
std::ostream&
operator<<(std::ostream& out, angle_map<T> const& m)
{
   for (typename angle_map<T>::const_iterator I = m.begin(); I != m.end(); ++I)
   {
      out << I->first << " " << I->second << '\n';
   }
   return out;
}

template <typename T>
angle_map<T>
operator+=(angle_map<T>& Input1, angle_map<T> const& Input2)
{
   for (auto const& I : Input2)
      Input1[I.first] += I.second;
   return Input1;
}

template <typename T>
angle_map<T>
operator-=(angle_map<T>& Input1, angle_map<T> const& Input2)
{
   for (auto const& I : Input2)
      Input1[I.first] -= I.second;
   return Input1;
}

template <typename T>
angle_map<T>
operator+(angle_map<T> const& Input1, angle_map<T> const& Input2)
{
   angle_map<T> Result(Input1);
   Result += Input2;
   return Result;
}

template <typename T>
angle_map<T>
operator-(angle_map<T> const& Input1, angle_map<T> const& Input2)
{
   angle_map<T> Result(Input1);
   Result -= Input2;
   return Result;
}

template <typename T, typename Scalar>
angle_map<T>
operator*=(angle_map<T>& Input, Scalar const& x)
{
   for (auto& I : Input)
      I.second *= x;
   return Input;
}

template <typename T, typename Scalar>
angle_map<T>
operator*(angle_map<T> Input, Scalar const& x)
{
   for (auto& I : Input)
      I.second = I.second * x;
   return Input;
}

template <typename T, typename Scalar>
angle_map<T>
operator*(Scalar const& x, angle_map<T> Input)
{
   for (auto& I : Input)
      I.second = x * I.second;
   return Input;
}

// Operators for vectors of angle maps.

template <typename T>
std::vector<angle_map<T>>
operator+=(std::vector<angle_map<T>>& Input1, std::vector<angle_map<T>> const& Input2)
{
   CHECK_EQUAL(Input1.size(), Input2.size());
   auto I1 = Input1.begin();
   auto I2 = Input2.begin();
   while (I1 != Input1.end())
   {
      *I1 += *I2;
      ++I1, ++I2;
   }
   return Input1;
}

template <typename T, typename Scalar>
std::vector<angle_map<T>>
ScalarMultiply(Scalar const& x, std::vector<angle_map<T>> Input)
{
   for (auto& I : Input)
      I = x * I;
   return Input;
}

#endif
