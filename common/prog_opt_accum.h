// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/prog_opt_accum.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
//
// accumulating_value for boost::program_options library
// Slightly modified from the version posted to the boost mailing list
// From: Bryan Green <bgreen@nas.nasa.gov>
// Date: Wed, 18 Jul 2007 15:28:59 -0700
// Message-ID: <200707182228.l6IMSxFv013280@ece06.nas.nasa.gov>
//
// Example:
//
// desc.add_options()
//    ("verbose,v", po_ext::accum_value<int>(), "print extra information")
//    ;
// ...
//    if (vm.count("verbose")) {
//        verbose = vm["verbose"].as<int>();
//
// It also supports the 'default_value()' method to initialize the counter.
//

#if !defined(PROG_OPT_ACCUM_H_SDJHCY43YH7O8YTOY)
#define PROG_OPT_ACCUM_H_SDJHCY43YH7O8YTOY

#include <boost/program_options.hpp>

namespace prog_opt_ext
{
   namespace prog_opt = boost::program_options;

   template <class T, class charT = char>
   class accumulating_value : public prog_opt::typed_value<T,charT>
   {
      public:
         accumulating_value(T* store_to=0)
         : prog_opt::typed_value<T,charT>(store_to), origin(0)
         {
            prog_opt::typed_value<T,charT>::zero_tokens();
         }

         accumulating_value* default_value(const T &v)
         {
            // setting a default value sets the origin to that value
            origin = v;
            prog_opt::typed_value<T,charT>::default_value(v);
            return this;
         }

         accumulating_value* default_value(const T &v,const std::string& textual)
         {
            // setting a default value sets the origin to that value
            origin = v;
            prog_opt::typed_value<T,charT>::default_value(v, textual);
            return this;
         }

         void xparse(boost::any& value_store,
                     const std::vector<std::basic_string<charT> >& new_tokens) const
         {
            // if this is the first occurrence of the option, initialize it
            // to the origin.
            if (value_store.empty())
               value_store = boost::any(origin);
            ++boost::any_cast<T&>(value_store);
         }

      private:
         T origin; // the numeric origin from which to increment upward.
   };

   template <class T>
   accumulating_value<T>*
   accum_value(T *v)
   {
      accumulating_value<T>* r =  new accumulating_value<T>(v);
      return r;
   }

   template <class T>
   accumulating_value<T>*
   accum_value()
   {
      return accum_value<T>(0);
   }

} // namespace prog_opt_ext

#endif
