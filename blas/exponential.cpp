// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/exponential.cpp
//
// Copyright (C) 2004-2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// A simple dense matrix class designed for scalar types.
//

#include "exponential.h"
#include "expokit/expokitf.h"

namespace blas
{

namespace detail
{
void Exponentiate(double t, int Size, std::complex<double> const* H, int ldH,
                  std::complex<double>* R, int ldR)
{
   // EXPOKIT doesn't handle matrices close to zero very well.  in this case, we
   // do a 1st order taylor series approximation
   double Norm = matrix_norm_frob_sq(Size, Size, H, ldH);
   if (Norm < 1E-10)
   {
      // 1st order approximation is R = 1+t*H
      matrix_copy_scaled('N', Size, Size, t, H, ldH, R, ldR);
      for (int i = 0; i < Size; ++i)
      {
	 R[(ldR+1)*i] += 1;
      }
      return;
   }

   // 2015-06-26: Increased degree of Pade approximation from 6 to 10, to get
   // better precision for low-dimensional operators
   // 2015-08-15: The sweet spot for accuracy of exponentials appears to be degree 9
   int ideg = 9;
   int lwork = 4*Size*Size + ideg + 1;
   std::complex<double>* work = new std::complex<double>[lwork];
   int* ipiv = new int[Size];
   int iexph;
   int ns;
   int iflag;

   EXPOKIT::zgpadm(ideg, Size, t, H, ldH, work, lwork, ipiv, iexph, ns, iflag);
   CHECK(iflag == 0)("EXPOKIT::zgpadm")(iflag);

   memcpy(R, work+iexph-1, Size*Size*sizeof(std::complex<double>));

   delete[] ipiv;
   delete[] work;
}

} // namespace detail

Matrix<std::complex<double>>
exp(Matrix<std::complex<double>> const& m)
{
   CHECK_EQUAL(m.rows(), m.cols());
   Matrix<std::complex<double>> Result(m.rows(), m.cols());

   detail::Exponentiate(1.0, m.rows(), m.storage(), m.leading_dimension(),
			Result.storage(), Result.leading_dimension());

   return Result;
}

} // namepsace blas
