// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/exponential.cpp
//
// Copyright (C) 2006-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
//
// implementation of eigen.h EXPOKIT wrappers
//
// Created 2006-05-12 Ian McCulloch
//

#include "eigen.h"
#include "expokit/expokitf.h"
#include "common/stackallocator.h"
#include <algorithm>

namespace LinearAlgebra
{

namespace Private
{

void Exponentiate(double t, int Size, std::complex<double> const* H, int ldH, 
                  std::complex<double>* R, int ldR)
{
   // 2015-06-26: Increased degree of Pade approximation from 6 to 10, to get
   // better precision for low-dimensional operators
   // 2015-08-15: The sweet spot for accuracy of exponentials appears to be degree 9
   int ideg = 9;
   int lwork = 4*Size*Size + ideg + 1;
   std::complex<double>* work = StackAlloc::allocate_type<std::complex<double> >(lwork);
   int* ipiv = StackAlloc::allocate_type<int>(Size);
   int iexph;
   int ns;
   int iflag;

   EXPOKIT::zgpadm(ideg, Size, t, H, ldH, work, lwork, ipiv, iexph, ns, iflag);
   CHECK(iflag == 0)("EXPOKIT::zgpadm")(iflag);

   memcpy(R, work+iexph-1, Size*Size*sizeof(std::complex<double>));

   StackAlloc::deallocate_type<int>(ipiv, Size);
   StackAlloc::deallocate_type<std::complex<double> >(work, lwork);
}

} // namespace Private

} // namespace LinearAlgebra
