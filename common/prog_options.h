// -*- C++ -*- $Id$
//
// Include wrapper for boost program_options library, adds FormatDefault() function
//

#if !defined(PROG_OPTIONS_H_SDCH348957Y23489J89A)
#define PROG_OPTIONS_H_SDCH348957Y23489J89A

#include "common/prog_opt_accum.h"
#include <string>
#include <boost/lexical_cast.hpp>

template <typename T>
std::string FormatDefault(std::string const& Text, T const& Value)
{
   return Text + " [default " + boost::lexical_cast<std::string>(Value) + "]";
}

#endif

