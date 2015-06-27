// -*- C++ -*- $Id$
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
   int ideg = 10;
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
