// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// interface/attributes.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(ATTRIBUTES_H_JKSHCUI87Y78HF784Y78YT7O8TO)
#define ATTRIBUTES_H_JKSHCUI87Y78HF784Y78YT7O8TO

#if defined(HAVE_CONFIG_H)  // we should include config.h before other headers
#include "config.h"
#endif
#include <map>
#include <string>
#include <ostream>
#include "pstream/pstream.h"

class Attribute
{
   public:
      Attribute();

      Attribute(Attribute const& a) : value_(a.value_) {}

      template <typename T>
      Attribute(T const& x);

      Attribute& operator=(Attribute const& a) { value_ = a.value_; return *this; }

      template <typename T>
      Attribute& operator=(T const& x);

      template <typename T>
      Attribute& operator+=(T const& x);

      bool empty() const { return value_.empty(); }

   //template <typename T>
   //operator T() const;

      template <typename T>
      T as() const;

      template <typename T>
      T get_or_default(T const& Default) const;

   private:
      std::string value_;

   friend std::ostream& operator<<(std::ostream& out, Attribute const& a);

   friend PStream::opstream& operator<<(PStream::opstream& out, Attribute const& a);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, Attribute& a);

};

class AttributeList
{
   private:
      typedef std::map<std::string, Attribute> data_type;

   public:
      typedef data_type::value_type value_type;
      typedef data_type::iterator iterator;
      typedef data_type::const_iterator const_iterator;

      AttributeList() {}

      AttributeList(AttributeList const& l) : data_(l.data_) {}

      Attribute operator[](std::string const& s) const;
      Attribute& operator[](std::string const& s) { return data_[s]; }

      int count(std::string const& s) const { return data_.count(s); }

      iterator begin() { return data_.begin(); }
      iterator end() { return data_.end(); }

      const_iterator begin() const { return data_.begin(); }
      const_iterator end() const { return data_.end(); }


   private:
      data_type data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, AttributeList const& l);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, AttributeList& l);
};

std::ostream& operator<<(std::ostream& out, AttributeList const& l);

#include "attributes.cc"

#endif
