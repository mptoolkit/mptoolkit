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

inline std::string ToStr(WhichEigenvalues w)
{
   switch (w)
   {
      case WhichEigenvalues::LargestMagnitude : return "LM";
      case WhichEigenvalues::SmallestMagnitude : return "SM";
      case WhichEigenvalues::LargestReal : return "LR";
      case WhichEigenvalues::SmallestReal : return "SR";
      case WhichEigenvalues::LargestImag : return "LI";
      case WhichEigenvalues::SmallestImag : return "SI";
      case WhichEigenvalues::LargestAlgebraic : return "LA";
      case WhichEigenvalues::SmallestAlgebraic : return "SA";
      case WhichEigenvalues::BothEnds : return "BE";
   }
   return "(unknown)";
}

struct CompareEigenvalues
{
   CompareEigenvalues(WhichEigenvalues w_) : w(w_) {}

   template <typename T>
   bool operator()(T a, T b)
   {
      switch (w)
      {
         case WhichEigenvalues::LargestMagnitude : return std::abs(a) > std::abs(b);
         case WhichEigenvalues::SmallestMagnitude : return std::abs(a) < std::abs(b);
         case WhichEigenvalues::LargestReal : return std::real(a) > std::real(b);
         case WhichEigenvalues::SmallestReal : return std::real(a) < std::real(b);
         case WhichEigenvalues::LargestImag : return false;
         case WhichEigenvalues::SmallestImag : return false;
         case WhichEigenvalues::LargestAlgebraic : return std::real(a) > std::real(b);
         case WhichEigenvalues::SmallestAlgebraic : return std::real(a) < std::real(b);
         case WhichEigenvalues::BothEnds : return a < b; // Do not sort.
      }
      return false;
   }

   template <typename T>
   bool operator()(std::complex<T> a, std::complex<T> b)
   {
      switch (w)
      {
         case WhichEigenvalues::LargestMagnitude : return std::abs(a) > std::abs(b);
         case WhichEigenvalues::SmallestMagnitude : return std::abs(a) < std::abs(b);
         case WhichEigenvalues::LargestReal : return std::real(a) > std::real(b);
         case WhichEigenvalues::SmallestReal : return std::real(a) < std::real(b);
         case WhichEigenvalues::LargestImag : return std::imag(a) > std::imag(b);
         case WhichEigenvalues::SmallestImag : return std::imag(a) < std::imag(b);
         case WhichEigenvalues::LargestAlgebraic : return std::real(a) > std::real(b);
         case WhichEigenvalues::SmallestAlgebraic : return std::real(a) < std::real(b);
         case WhichEigenvalues::BothEnds : PANIC("undefined"); return false; // Do not sort.
      }
      PANIC("invalid");
      return false;
   }

   WhichEigenvalues w;
};

template <typename MultFunc>
Vector<std::complex<double>>
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, WhichEigenvalues which, std::complex<double> const* InitialGuess, double tol,
                  std::vector<std::complex<double>>* OutputVectors,
                  int ncv, bool Sort, int Verbose)
{
   std::set<WhichEigenvalues> ValidEigenvaluesComplex = { WhichEigenvalues::LargestMagnitude, WhichEigenvalues::SmallestMagnitude, WhichEigenvalues::LargestReal, WhichEigenvalues::SmallestReal, WhichEigenvalues::LargestImag, WhichEigenvalues::SmallestImag };

   CHECK(ValidEigenvaluesComplex.count(which))("Invalid eigenvalue selection for DiagonalizeARPACK<complex>()");
   if (Verbose >= 1)
   {
      std::cerr << "Total dimension = " << n << std::endl;
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
      // 2022-02-17: this bug isn't actually harmless - it means that the caller can get
      // more elements in the array than expected.  Truncate the eigenvalue array if we get more than expected.
      if (Verbose >= 1)
      {
         std::cerr << "Constructing matrix for direct diagonalization" << std::endl;
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
         (*OutputVectors) = std::vector<std::complex<double>>(n*n);
         for (int k =0; k < n; ++k)
         {
            LinearAlgebra::make_vec(&(*OutputVectors)[n*k], n) = RV(k, LinearAlgebra::all);
         }
      }
      // If we are using LAPACK for the diagonalization then we might have more eigenvalues than we need.  In order to ensure
      // that we get the right eigenvalues, we must sort them, even if the caller doesn't require it.
      Sort = true;
   }
   else
   {
      // arpack parameters
      int ido = 0;  // first call
      char bmat = 'I'; // standard eigenvalue problem
      std::string w = ToStr(which);
      int const nev = std::min(NumEigen, n-2); // number of eigenvalues to be computed
      std::vector<std::complex<double> > resid(n);  // residual
      ncv = std::min(std::max(ncv, 2*nev + 10), n);            // length of the arnoldi sequence
      std::vector<std::complex<double> > v(n*ncv);   // arnoldi vectors
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
      if (InitialGuess)
      {
         std::copy(InitialGuess, InitialGuess+n, v.data());
      }
      int info = InitialGuess ? 1 : 0;  // this is the indicator whether to use the initial residual, or create one randomly

      if (Verbose >= 1)
      {
         std::cerr << "Starting ARPACK mode " << w << std::endl;
         if (Verbose >= 2)
            std::cerr << "n=" << n << ", nev=" << nev << ", ncv=" << ncv << std::endl;
      }

      int NumMultiplies = 0;

      ARPACK::znaupd(&ido, bmat, n, w.c_str(), nev, tol, &resid[0], ncv,
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

         ARPACK::znaupd(&ido, bmat, n, w.c_str(), nev, tol, &resid[0], ncv,
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
         std::cerr << "Finished ARPACK, nev=" << nev << ", ncv=" << ncv << ", NumMultiplies=" << NumMultiplies << " " << iparam << std::endl;
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
                     bmat, n, w.c_str(), nev, tol, &resid[0], ncv, &v[0], ldv,
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
      CompareEigenvalues C(which);
      // a simple exchange sort
      for (unsigned i = 0; i < Result.size()-1; ++i)
      {
         for (unsigned j = i+1; j < Result.size(); ++j)
         {
            if (C(Result[j], Result[i]))
            {
               std::swap(Result[i], Result[j]);
               if (OutputVectors)
               {
                  std::vector<std::complex<double>> Temp(&(*OutputVectors)[n*i],
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

   // If we ended up with more eigenvalues than we wanted, resize the arrays.
   if (Result.size() > NumEigen)
   {
      LinearAlgebra::Vector<std::complex<double>> R2 = Result[LinearAlgebra::range(0,NumEigen)];
      std::swap(Result, R2);
      if (OutputVectors)
         OutputVectors->resize(n*NumEigen);
   }

   // restore the ARPACK debug log level before returning
   ARPACK::debug().mceupd = DebugOutputLevelSave;

   return Result;
}

} // namespace LinearAlgebra
