// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/spectrum_arpack.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Get the spectrum of an operator using ARPACK

#include "spectrum_arpack.h"


struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_)
      : L(L_), QShift(QShift_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
         r = operator_prod(herm(*I), r, *I);
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
};

LinearAlgebra::Vector<std::complex<double> >
get_spectrum(LinearWavefunction const& Psi, QuantumNumber const& QShift, int NumEigen,
             QuantumNumbers::QuantumNumber const& q, double tol,
             LinearAlgebra::Vector<MatrixOperator>* OutputVectors, int ncv, int Verbose)
{
   PackMatrixOperator Pack(Psi.Basis2(), Psi.Basis2(), q);

   // arpack parameters
   int ido = 0;  // first call
   char bmat = 'I'; // standard eigenvalue problem
   int const n = Pack.size();
   char which[3] = "LR";                      // largest real part
   //   char which[3] = "LM";                      // largest magnitude
   int const nev = std::min(NumEigen, n-2);   // number of eigenvalues to be computed
   std::vector<std::complex<double> > resid(n);  // residual
   ncv = std::min(std::max(ncv, 2*nev), n);            // length of the arnoldi sequence
   std::vector<std::complex<double> > v(n*ncv);   // arnoldi vectors
   int const ldv = n;
   ARPACK::iparam_t iparam;
   iparam.ishift = 1;      // exact shifts
   iparam.mxiter = 10000;  // maximum number of arnoldi iterations (restarts?)
   iparam.mode = 1;  // ordinery eigenvalue problem
   ARPACK::zn_ipntr_t ipntr;
   std::vector<std::complex<double> > workd(3*n);
   int const lworkl = 3*ncv*ncv + 5*ncv;
   std::vector<std::complex<double> > workl(lworkl);
   std::vector<double> rwork(ncv);
   int info = 0;  // no initial residual

   LeftMultiply Mult(Psi, QShift);  // multiply functor

   if (Verbose >= 1)
   {
      std::cerr << "Total dimension = " << n << '\n';
      std::cerr << "Starting ARPACK";
   }

   int NumMultiplies = 0;

   //   TRACE(tol);

   ARPACK::znaupd(&ido, bmat, n, which, nev, tol, &resid[0], ncv,
                  &v[0], ldv, &iparam, &ipntr, &workd[0],
                  &workl[0], lworkl, &rwork[0], &info);
   CHECK(info >= 0)(info);

   while (ido != 99)
   {
      if (ido == -1 || ido == 1)
      {
         MatrixOperator x = Pack.unpack(&workd[ipntr.x]);
         x = Mult(x);
         if (Verbose >= 2)
            std::cerr << '.';
         Pack.pack(x, &workd[ipntr.y]);
         ++NumMultiplies;
      }
      else
      {
         PANIC("unexpected reverse communication operation.")(ido);
      }

      ARPACK::znaupd(&ido, bmat, n, which, nev, tol, &resid[0], ncv,
                     &v[0], ldv, &iparam, &ipntr, &workd[0],
                     &workl[0], lworkl, &rwork[0], &info);
      CHECK(info >= 0)(info);
   }

   if (Verbose >= 1)
   {
      std::cerr << "\nFinished ARPACK, nev=" << nev << ", ncv=" << ncv << ", NumMultiplies=" << NumMultiplies << '\n';
   }

   // get the eigenvalues
   bool rvec = OutputVectors != NULL; // should we calculate eigenvectors?
   char howmny = 'A';
   std::vector<int> select(ncv);
   std::vector<std::complex<double> > d(nev+1);
   std::vector<std::complex<double> > z(OutputVectors ? n*(nev+1) : 1); // output array
   int ldz = n;
   std::complex<double> sigma;   // not referenced
   std::vector<std::complex<double> > workev(2*ncv);
   ARPACK::zneupd(rvec, howmny, &select[0], &d[0], &z[0], ldz, sigma, &workev[0],
                  bmat, n, which, nev, tol, &resid[0], ncv, &v[0], ldv,
                  &iparam, &ipntr, &workd[0],
                  &workl[0], lworkl, &rwork[0], &info);
   CHECK(info >= 0)(info);

   LinearAlgebra::Vector<std::complex<double> > Result(nev);
   for (int i = 0; i < nev; ++i)
   {
      Result[i] = d[i];
   }

   // eigenvectors
   if (OutputVectors)
   {
      *OutputVectors = LinearAlgebra::Vector<MatrixOperator>(nev);
      for (int i = 0; i < nev; ++i)
      {
         (*OutputVectors)[i] = Pack.unpack(&z[n*i]);
      }
   }

   return Result;
}
