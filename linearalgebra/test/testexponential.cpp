// -*- C++ -*- $Id$
//
// test the expokit exponential function.
// The main thing we want to do is see if it works for small inputs.

#include "linearalgebra/exponential.h"
#include "linearalgebra/scalarmatrix.h"


using namespace LinearAlgebra;

int main()
{
   int const sz = 100;
   Matrix<std::complex<double> > M(sz,sz,0.0);

   M = ScalarMatrix<double>(sz, sz, 1.0);

   double x = 1.0;
   double e = exp(1/x);
   double r = norm_frob(Exponentiate(1.0, Matrix<std::complex<double> >((1/x) * M)));
   while (norm_frob(r/std::sqrt(double(sz)) - e) < 1e-10)
   {
      x = x * 2;
      e = exp(1/x);
      r = norm_frob(Exponentiate(1.0, Matrix<std::complex<double> >((1/x) * M)));
      TRACE(x)(e);
   }
   TRACE(x)(e)(r);
}
