// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// expokit/expokitf.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  expokitf.h

  C++ interface to EXPOKIT.
*/

#if !defined(EXPOKITF_H_SDFHJK34893R489UE8ER89FDJOIJ)
#define EXPOKITF_H_SDFHJK34893R489UE8ER89FDJOIJ

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include "common/fortran.h"
#include "common/trace.h"
#include <complex>

#if defined(EXPOKIT_TRACE_DETAILED)
#define TRACE_EXPOKIT(Msg) TRACE(Msg)
#else
#define TRACE_EXPOKIT(Msg) DUMMY_TRACE(Msg)
#endif

namespace EXPOKIT
{

using namespace Fortran;

// wrapper functions

// complex

void zgpadm(integer ideg, integer m, double t, std::complex<double> const* H, integer ldh,
            std::complex<double>* restrict work, integer lwork, integer* restrict ipiv,
            integer& restrict iexph, integer& restrict ns, integer& restrict iflag);

namespace raw
{

extern "C"
{

void F77NAME(zgpadm)(integer const* ideg, integer const* m, double const* t,
                     complex const* H, integer const* ldh,
                     complex* restrict work, integer const* lwork,
                     integer* restrict ipiv, integer* restrict iexph,
                     integer* restrict ns, integer* restrict iflag);

} // extern "C"

} // namespace raw

// implementation of the wrappers

void zgpadm(integer ideg, integer m, double t, std::complex<double> const* H, integer ldh,
            std::complex<double>* restrict work, integer lwork, integer* restrict ipiv,
            integer& restrict iexph, integer& restrict ns, integer& restrict iflag)
{
   TRACE_EXPOKIT("zgpadm")(ideg)(m)(t)(H)(ldh)(work)(lwork)(ipiv)(iexph)(ns)(iflag);
   raw::F77NAME(zgpadm)(&ideg, &m, &t, reinterpret_cast<complex const*>(H),
                        &ldh, reinterpret_cast<complex*>(work), &lwork,
                        ipiv, &iexph, &ns, &iflag);
}

} // namespace EXPOKIT

/*
      subroutine ZGPADM(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)

      implicit none
      double precision t
      integer          ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
      complex*16       H(ldh,m), wsp(lwsp)

*-----Purpose----------------------------------------------------------|
*
*     Computes exp(t*H), the matrix exponential of a general complex
*     matrix in full, using the irreducible rational Pade approximation
*     to the exponential exp(z) = r(z) = (+/-)( I + 2*(q(z)/p(z)) ),
*     combined with scaling-and-squaring.
*
*-----Arguments--------------------------------------------------------|
*
*     ideg      : (input) the degre of the diagonal Pade to be used.
*                 a value of 6 is generally satisfactory.
*
*     m         : (input) order of H.
*
*     H(ldh,m)  : (input) argument matrix.
*
*     t         : (input) time-scale (can be < 0).
*
*     wsp(lwsp) : (workspace/output) lwsp .ge. 4*m*m+ideg+1.
*
*     ipiv(m)   : (workspace)
*
*>>>> iexph     : (output) number such that wsp(iexph) points to exp(tH)
*                 i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
*                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
*                 NOTE: if the routine was called with wsp(iptr),
*                       then exp(tH) will start at wsp(iptr+iexph-1).
*
*     ns        : (output) number of scaling-squaring used.
*
*     iflag     : (output) exit flag.
*                       0 - no problem
*                      <0 - problem
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|
*
*/

#endif
