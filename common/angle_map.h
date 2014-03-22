// -*- C++ -*- $Id$

// angle_map.h
// a std::map for angular quantities

#if !defined(MPTOOLKIT_COMMON_ANGLE_MAP_H)
#define MPTOOLKIT_COMMON_ANGLE_MAP_H

#include <complex>
#include <map>
#include <cmath>

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
      TRACE("Adding new angle")(c);
      return Map[c];
   }
   // else
   return I->second;
}





#endif
