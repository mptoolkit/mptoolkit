// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
#include "common/math_const.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> v = 0.5 * math_const::pi * Range(0, 10);
   Vector<double> w = sin(v);
   CHECK_CLOSE(w[0], 0);
   CHECK_CLOSE(w[1], 1);
   CHECK_CLOSE(w[2], 0);
   CHECK_CLOSE(w[3], -1);
   CHECK_CLOSE(norm_1(w), 5);

   w = abs(w);
   CHECK_CLOSE(w[0], 0);
   CHECK_CLOSE(w[1], 1);
   CHECK_CLOSE(w[2], 0);
   CHECK_CLOSE(w[3], 1);
   CHECK_CLOSE(norm_1(w), 5);

   w = cos(0.5 * math_const::pi * Range(0, 10));
   CHECK_CLOSE(w[0], 1);
   CHECK_CLOSE(w[1], 0);
   CHECK_CLOSE(w[2], -1);
   CHECK_CLOSE(w[3], 0);
   CHECK_CLOSE(norm_1(w), 5);
   
   w = exp(w);
   CHECK_CLOSE(w[0], std::exp(1.0));
   CHECK_CLOSE(w[1], 1);
   CHECK_CLOSE(w[2], std::exp(-1.0));
   CHECK_CLOSE(w[3], 1);

   Vector<std::complex<double> > cv(2);
   cv[0] = std::complex<double>(1.0, 0.0);
   cv[1] = std::complex<double>(0.0, 2.0);

   Vector<double> cva = abs(cv);
   CHECK_CLOSE(cva[0], 1);
   CHECK_CLOSE(cva[1], 2);

   Vector<double> cvarg = arg(cv);
   CHECK_CLOSE(cvarg[0], 0);
   CHECK_CLOSE(cvarg[1], 0.5 * math_const::pi);

   Vector<std::complex<double> > xv = std::complex<double>(0.0,0.5) * (math_const::pi * Range(0, 10));
   cv = exp(xv);
   CHECK_CLOSE(cv[0], 1.0);
   CHECK_CLOSE(cv[1], std::complex<double>(0,1));
   CHECK_CLOSE(cv[2], -1.0);
   CHECK_CLOSE(cv[3], std::complex<double>(0,-1));
   CHECK_CLOSE(norm_1(cv), 10);
}
