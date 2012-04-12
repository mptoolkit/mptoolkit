// -*- C++ -*- $Id$

#include <boost/utility/result_of.hpp>

#include <iostream>
#include <iomanip>

#include "common/trace.h"

using tracer::typeid_name;

int F(int, int) { return 0; };

struct f
{
   int operator()(int, int) const { return 1; };

   template <typename T> struct result;

   template <typename T, typename U>
   struct result<f(T const&,U const&)> { typedef float type; };
};

int main()
{
   TRACE(typeid_name(boost::result_of<f(int,int)>::type()));
}
