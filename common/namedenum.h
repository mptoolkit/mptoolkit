// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/namedenum.h
//
// Copyright (C) 2024 Ian McCulloch <ian@qusim.net>
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

// A class to represent an enumeration that is iterable and has a string name
// associated with each item.
//
// To use, define a class or struct that has 4 members:
// Enum        - must be an enumeration type
// Default     - a constexpr value of type Enum, that is the value of a default-constructed NamedEnumeration
// StaticName  - must be a static constexpr of type char const* that gives a desciption of the enumeration.
// Names       - must be a static constexpr array of strings, of exactly the same size as Enum.
//
// An example:
// struct MyEnumTraits
// {
//    enum Enum { Some, Enumeration, Elements };
//    static constexpr Enum Default = Enumeration;
//    static constexpr char const* StaticName = "the example enumeration";
//    static constexpr std::array<char const*, 3> Names = { "some", "enumeration", "elements" };
// };
//

#include <array>
#include <string>
#include <exception>

#if !defined(MPTOOLKIT_COMMON_NAMEDENUM_H)
#define MPTOOLKIT_COMMON_NAMEDENUM_H

// An enumeration that supports iteration, as well as
template <typename Traits>
class NamedEnumeration : public Traits
{
   public:
      using Enum = typename Traits::Enum;
      using Traits::StaticName;
      static constexpr std::size_t N = Traits::Names.size();
      static constexpr Enum DEFAULT = Traits::Default;
      static constexpr Enum BEGIN = static_cast<Enum>(0);
      static constexpr Enum END = static_cast<Enum>(N);

      NamedEnumeration() : e(DEFAULT) {}

      NamedEnumeration(Enum a) : e(a) {}

      explicit NamedEnumeration(std::string Name);

      // Enable iteration (including range-based for loop) over the available algorithms
      NamedEnumeration begin() const { return BEGIN; }
      NamedEnumeration end() const { return END; }

      static constexpr std::size_t size() { return N; }

      bool operator==(NamedEnumeration const& Other) const { return e == Other.e; }
      bool operator!=(NamedEnumeration const& Other) const { return e != Other.e; }
      bool operator==(Enum a) const { return e == a; }
      bool operator!=(Enum a) const { return e != a; }
      NamedEnumeration& operator++() { e = static_cast<Enum>(e+1); return *this; }
      NamedEnumeration& operator--() { e = static_cast<Enum>(e-1); return *this; }
      const NamedEnumeration& operator*() const { return *this; }

      static std::string ListAvailable();

      std::string Name() const { return Traits::Names[e]; }

   private:
      Enum e;
};

template <typename Traits>
std::string NamedEnumeration<Traits>::ListAvailable()
{
   std::string Result;
   bool first = true;
   for (auto a : NamedEnumeration())
   {
      if (!first)
         Result += ", ";
      Result += a.Name();
      first = false;
   }
   return Result;
}

template <typename Traits>
NamedEnumeration<Traits>::NamedEnumeration(std::string Name)
{
   std::transform(Name.begin(), Name.end(), Name.begin(), [](unsigned char c){ return std::tolower(c); });
   for (auto a : NamedEnumeration())
   {
      if (a.Name() == Name)
      {
         e = a.e;
         return;
      }
   }
   using namespace std::literals::string_literals;
   std::string ErrorStr = "Unknown initializer for "s + StaticName + "; choices are " + this->ListAvailable() + '.';
   throw std::runtime_error(ErrorStr);
}

#endif
