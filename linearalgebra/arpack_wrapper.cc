// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/arpack_wrapper.cc
//
// Copyright (C) 2007-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "common/arpackf.h"
#include "linearalgebra/vectormemproxy.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"

namespace LinearAlgebra
{

template <typename MultFunc>
Vector<std::complex<double> >
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, double tol,
                  std::vector<std::complex<double>>* OutputVectors,
                  int ncv, bool Sort, int Verbose)
{
   if (Verbose >= 1)
   {
      std::cerr << "Total dimension = " << n << '\n';
   }
   // save the arpack debug log level so we can restore it later
   int DebugOutputLevelSave = ARPACK::debug().mceupd;
   ARPACK::debug().mceupd = Verbose;

   LinearAlgebra::Vector<std::complex<double>> Result;

   if (NumEigen >= n-1)
   {
      // This is a case that ARPACK doesn't handle - it is apparently not
      // capable of producting more than n-2 eigenvalues.  Instead we handle this
      // case by converting to a dense matrix.
      // 2014-04-01: n-2 because we require NCV-NEV >= 2, and NCV <= N.
      // There is a (harmless?) bug here in that if we only want n-1 eigenvalues
      // then we'll actually get all n of them.
      if (Verbose >= 1)
      {
         std::cerr << "Constructing matrix for direct diagonalization\n";
      }
      LinearAlgebra::Vector<std::complex<double>> e(n, 0.0), Out(n);
      LinearAlgebra::Matrix<std::complex<double>> Mat(n, n);
      for (int k = 0; k < n; ++k)
      {
         e[k] = 1.0;
         Mult(data(e), data(Out));
         Mat(LinearAlgebra::all, k) = Out;
         e[k] = 0.0;
      }
      LinearAlgebra::Matrix<std::complex<double>> LV, RV;
      Result = LinearAlgebra::Diagonalize(Mat, LV, RV);
      if (OutputVectors)
      {
         (*OutputVectors) = std::vector<std::complex<double> >(n*n);
         for (int k =0; k < n; ++k)
         {
            LinearAlgebra::make_vec(&(*OutputVectors)[n*k], n) = RV(k, LinearAlgebra::all);
         }
      }
   }
   else
   {
      // arpack parameters
      int ido = 0;  // first call
      char bmat = 'I'; // standard eigenvalue problem
      char which[3] = "LM";                      // largest magnitude
      int const nev = std::min(NumEigen, n-2); // number of eigenvalues to be computed
      std::vector<std::complex<double>> resid(n);  // residual
      ncv = std::min(std::max(ncv, 2*nev + 10), n);            // length of the arnoldi sequence
      std::vector<std::complex<double>> v(n*ncv);   // arnoldi vectors
      int const ldv = n;
      ARPACK::iparam_t iparam;
      iparam.ishift = 1;      // exact shifts
      iparam.mxiter = 10000;  // maximum number of arnoldi iterations (restarts?)
      iparam.mode = 1;  // ordinary eigenvalue problem
      ARPACK::zn_ipntr_t ipntr;
      std::vector<std::complex<double>> workd(3*n);
      int const lworkl = 3*ncv*ncv + 5*ncv;
      std::vector<std::complex<double>> workl(lworkl);
      std::vector<double> rwork(ncv);
      int info = 0;  // no initial residual

      if (Verbose >= 1)
      {
         std::cerr << "Starting ARPACK mode LM\n";
         if (Verbose >= 2)
            std::cerr << "n=" << n << ", nev=" << nev << ", ncv=" << ncv << '\n';
      }

      int NumMultiplies = 0;

      ARPACK::znaupd(&ido, bmat, n, which, nev, tol, &resid[0], ncv,
                     &v[0], ldv, &iparam, &ipntr, &workd[0],
                     &workl[0], lworkl, &rwork[0], &info);
      CHECK(info >= 0)(info)(n)(nev)(ncv);

      while (ido != 99)
      {
         if (ido == -1 || ido == 1)
         {
            if (Verbose >= 2)
               std::cerr << '.';
            Mult(&workd[ipntr.x], &workd[ipntr.y]);
            ++NumMultiplies;
         }
         else
         {
            PANIC("unexpected reverse communication operation.")(ido);
         }

         ARPACK::znaupd(&ido, bmat, n, which, nev, tol, &resid[0], ncv,
                        &v[0], ldv, &iparam, &ipntr, &workd[0],
                        &workl[0], lworkl, &rwork[0], &info);
         if (info == -9)
         {
            // info == -9 indicates that the starting vector is zero.  Since we have already done a matrix-vector
            // multiply, the only way this can happen is if the matrix itself is zero.  Hence we know what the
            // eigenvalues are.
            Result = LinearAlgebra::Vector<std::complex<double>>(nev, std::complex<double>(0.0, 0.0));
            if (OutputVectors)
            {
               *OutputVectors = std::vector<std::complex<double>>(n*nev);
               // construct an eigenbasis; this is arbitrary but we can choose the standard basis for the first
               // nev eigenvectors
               for (int i = 0; i < nev; ++i)
               {
                  (*OutputVectors)[i*n + i] = 1.0;
               }
            }
            return Result;
         }
         CHECK(info >= 0)(info)(n)(nev)(ncv);
      }

      if (Verbose >= 1)
      {
         std::cerr << "\nFinished ARPACK, nev=" << nev << ", ncv=" << ncv << ", NumMultiplies=" << NumMultiplies << " " << iparam << '\n';
      }

      // get the eigenvalues
      bool rvec = OutputVectors != nullptr; // should we calculate eigenvectors?
      char howmny = 'A';
      std::vector<int> select(ncv);
      std::vector<std::complex<double>> d(nev+1);
      std::vector<std::complex<double>> z(OutputVectors ? n*nev : 1); // output array
      int ldz = n;
      std::complex<double> sigma;   // not referenced
      std::vector<std::complex<double>> workev(2*ncv);
      ARPACK::zneupd(rvec, howmny, &select[0], &d[0], &z[0], ldz, sigma, &workev[0],
                     bmat, n, which, nev, tol, &resid[0], ncv, &v[0], ldv,
                     &iparam, &ipntr, &workd[0],
                     &workl[0], lworkl, &rwork[0], &info);
      CHECK(info >= 0)("arpack::zneupd")(info)(nev)(ncv);

      Result = LinearAlgebra::Vector<std::complex<double>>(nev);
      for (int i = 0; i < nev; ++i)
      {
         Result[i] = d[i];
      }

      // eigenvectors
      if (OutputVectors)
      {
         OutputVectors->empty();
         std::swap(z, *OutputVectors);
      }

   }

   if (Sort)
   {
      // a simple exchange sort
      for (unsigned i = 0; i < Result.size()-1; ++i)
      {
         for (unsigned j = i+1; j < Result.size(); ++j)
         {
            if (norm_frob(Result[j]) > norm_frob(Result[i]))
            {
               std::swap(Result[i], Result[j]);
               if (OutputVectors)
               {
                  std::vector<std::complex<double> > Temp(&(*OutputVectors)[n*i],
                                                          &(*OutputVectors)[n*i]+n);
                  LinearAlgebra::fast_copy(&(*OutputVectors)[n*j],
                                           &(*OutputVectors)[n*j]+n,
                                           &(*OutputVectors)[n*i]);
                  LinearAlgebra::fast_copy(&Temp[0], &Temp[0]+n, &(*OutputVectors)[n*j]);
               }
            }
         }
      }
   }

   // restore the ARPACK debug log level before returning
   ARPACK::debug().mceupd = DebugOutputLevelSave;

   return Result;
}

} // namespace LinearAlgebra
