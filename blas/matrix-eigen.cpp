// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/eigen.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// implementation of matrix-eigen.h LAPACK wrappers
//
// Created 2004-06-07 Ian McCulloch
//

#include "matrix-eigen.h"
#include "common/lapackf.h"

namespace blas
{

namespace detail
{

void DiagonalizeSymmetric(int Size, double* Data, int LeadingDim, double* Eigen)
{
   char jobz = 'V';
   char uplo = 'L';
   double worksize;
   double* work = &worksize;
   int lwork = -1;                         // for query of lwork
   Fortran::integer info = 0;

   // query for the optimial size of the workspace
   LAPACK::dsyev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, info);

   lwork = int(work[0]);
   work = new double[lwork];

   // do the actual call
   LAPACK::dsyev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, info);

   CHECK(info == 0)("LAPACK::dsyev")(info);

   delete[] work;
}

void DiagonalizeHermitian(int Size, std::complex<double>* Data, int LeadingDim, double* Eigen)
{
   if (Size == 0)
   {
      DEBUG_WARNING("Zero size DiagonalizeHermitian");
      return;
   }

   char jobz = 'V';
   char uplo = 'L';
   std::complex<double>  worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;                         // for query of lwork
   double* rwork;
   int lrwork = std::max(1, 3*Size-2);
   Fortran::integer info = 0;

   rwork = new double[lrwork];

   // query for the optimial size of the workspace
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   CHECK(info == 0)("LAPACK::zheev")(info);

   delete[] work;
   delete[] rwork;
}

} // namespace detail

} // namespace blas
