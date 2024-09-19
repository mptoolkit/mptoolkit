// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/blas1f.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ian@qusim.net>
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
/* -*- C++ -*-
  blas1f.h

  C++ interface to BLAS level 1

  Created 2000-09-22 Ian McCulloch

  Ideas from LAPACK++
*/

#if !defined(MPTOOLKIT_COMMON_BLAS1_H)
#define MPTOOLKIT_COMMON_BLAS1_H

#if defined(HAVE_CONFIG_H)
#include "config.h"   // to get restrict working
#else
#if !defined(restrict)
#define restrict
#endif
#endif

#include "fortran.h"
#include "trace.h"

#if defined(BLAS1_TRACE_DETAILED)
#define TRACE_BLAS1(Msg) TRACE(Msg)
#else
#define TRACE_BLAS1(Msg) DUMMY_TRACE(Msg)
#endif

namespace BLAS
{

using namespace Fortran;

// wrapper functions for BLAS level 1

// real

double dasum(integer n, double const* dx, integer incx);

void daxpy(integer n, double da, double const* dx, integer incx, double* restrict dy, integer incy);

void dcopy(integer n, double const* dx, integer incx, double* restrict dy, integer incy);

double ddot(integer n, double const* dx, integer incx, double const* dy, integer incy);

double dnrm2(integer n, double const* dx, integer incx);

void drot(integer n, double* restrict dx, integer incx,
          double* restrict dy, integer incy, double c, double s);

void drotg(double* restrict da, double* restrict db, double* restrict c, double* restrict s);

void dscal(integer n, double da, double* dx, integer incx);

void dswap(integer n, double* restrict dx, integer incx, double* restrict dy, integer incy);

// return value here is adjusted to be zero-based
integer idamax(integer n, double const* dx, integer incx);

// complex

double dzasum(integer n, std::complex<double> const* dx, integer incx);

void zaxpy(integer n, std::complex<double> da,
           std::complex<double> const* dx, integer incx,
           std::complex<double>* restrict dy, integer incy);

void zcopy(integer n, std::complex<double> const* dx, integer incx,
           std::complex<double>* restrict dy, integer incy);

#if defined(HAVE_FORTRAN_COMPLEX_RETURN)

std::complex<double> zdotu(integer n, std::complex<double> const* dx, integer incx,
                           std::complex<double> const* dy, integer incy);
#define HAVE_ZDOTU

std::complex<double> zdotc(integer n, std::complex<double> const* dx, integer incx,
                           std::complex<double> const* dy, integer incy);
#define HAVE_ZDOTC

#endif

double dznrm2(integer n, std::complex<double> const* dx, integer incx);

void zrot(integer n, std::complex<double>* restrict dx, integer incx,
          std::complex<double>* restrict dy, integer incy,
          std::complex<double> c, std::complex<double> s);

void zrotg(std::complex<double>* restrict da, std::complex<double>* restrict db,
           std::complex<double>* restrict c, std::complex<double>* restrict s);

void zscal(integer n, std::complex<double> da, std::complex<double>* dx, integer incx);

void zdscal(integer n, double da, std::complex<double>* dx, integer incx);

void zswap(integer n, std::complex<double>* restrict dx, integer incx,
           std::complex<double>* restrict dy, integer incy);

// return value here is adjusted to be zero-based
integer izamax(integer n, std::complex<double> const* dx, integer incx);

// the raw functions are in their own namespace, the wrapper functions call these.
namespace raw
{

extern "C"
{

// real

double F77NAME(dasum)(const integer *n, const double *dx, const integer *incx);

void F77NAME(daxpy)(const integer *n, const double *da, const double *dx,
                    const integer *incx, double* restrict dy, const integer *incy);

void F77NAME(dcopy)(const integer *n, double const* dx, const integer *incx, double *dy,
                    const integer *incy);

double F77NAME(ddot)(const integer *n, const double *dx, const integer *incx,
                     const double *dy, const integer *incy);

double F77NAME(dnrm2)(const integer *n, const double *dx, const integer *incx);

void F77NAME(drot)(const integer *n, double* restrict dx, integer const* incx, double* restrict dy,
                        integer const* incy, const double *c, const double *s);

void F77NAME(drotg)(double* restrict da, double* restrict db, double* restrict c, double* restrict s);

void F77NAME(dscal)(const integer *n, double* const da, double* dx, const integer *incx);

void F77NAME(dswap)(const integer *n, double *dx, const integer *incx, double *dy,
                        const integer *incy);

integer F77NAME(idamax)(const integer *n, const double *dx, const integer *incx);

// complex

double F77NAME(dzasum)(const integer *n, const std::complex<double> *dx, const integer *incx);

void F77NAME(zaxpy)(const integer *n, const std::complex<double> *da, const std::complex<double> *dx,
                    const integer *incx, std::complex<double>* restrict dy, const integer *incy);

void F77NAME(zcopy)(const integer *n,
                    std::complex<double> const* dx, const integer *incx,
                    std::complex<double> *dy, const integer *incy);

#if defined(FORTRAN_COMPLEX_RETURN_FIRST_ARG)
void F77NAME(zdotu)(std::complex<double>* result, const integer *n,
                    const std::complex<double> *dx, const integer *incx,
                    const std::complex<double> *dy, const integer *incy);

void F77NAME(zdotc)(std::complex<double>* result, const integer *n,
                    const std::complex<double> *dx, const integer *incx,
                    const std::complex<double> *dy, const integer *incy);
#elif defined(FORTRAN_COMPLEX_RETURN_IN_REGISTER)
Fortran::complex
F77NAME(zdotu)(const integer *n,
               const std::complex<double> *dx, const integer *incx,
               const std::complex<double> *dy, const integer *incy);
Fortran::complex
F77NAME(zdotc)(const integer *n,
               const std::complex<double> *dx, const integer *incx,
               const std::complex<double> *dy, const integer *incy);
#endif

double F77NAME(dznrm2)(const integer *n,
                       const std::complex<double> *dx, const integer *incx);

void F77NAME(zrot)(const integer *n, std::complex<double>* restrict dx, integer const* incx,
                   std::complex<double>* restrict dy, integer const* incy,
                   const std::complex<double> *c, const std::complex<double> *s);

void F77NAME(zrotg)(std::complex<double>* restrict da, std::complex<double>* restrict db,
                    std::complex<double>* restrict c, std::complex<double>* restrict s);

void F77NAME(zscal)(const integer *n, std::complex<double>* const da,
                    std::complex<double>* dx, const integer *incx);

void F77NAME(zdscal)(const integer *n, double* const da,
                     std::complex<double>* dx, const integer *incx);

void F77NAME(zswap)(const integer *n, std::complex<double> *dx, const integer *incx,
                    std::complex<double> *dy, const integer *incy);

integer F77NAME(izamax)(const integer *n, const std::complex<double> *dx, const integer *incx);

} // extern "C"

} // namespace raw

// inlines

// real

inline
double dasum(integer n, double const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: dasum")(n)(dx)(incx);
   return raw::F77NAME(dasum)(&n, dx, &incx);
}

inline
void daxpy(integer n, double da, double const* dx, integer incx, double* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: daxpy")(n)(da)(dx)(incx)(dy)(incy);
   raw::F77NAME(daxpy)(&n, &da, dx, &incx, dy, &incy);
}

inline
void dcopy(integer n, double const* dx, integer incx, double* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: dcopy")(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(dcopy)(&n, dx, &incx, dy, &incy);
}

inline
double ddot(integer n, double const* dx, integer incx, double const* dy, integer incy)
{
   TRACE_BLAS1("BLAS1: ddot")(n)(dx)(incx)(dy)(incy);
   return raw::F77NAME(ddot)(&n, dx, &incx, dy, &incy);
}

inline
double dnrm2(integer n, double const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: dnrm2")(n)(dx)(incx);
   return raw::F77NAME(dnrm2)(&n, dx, &incx);
}

inline
void drot(integer n, double* restrict dx, integer incx,
          double* restrict dy, integer incy, double c, double s)
{
   TRACE_BLAS1("BLAS1: drot")(n)(dx)(incx)(dy)(incy)(c)(s);
   raw::F77NAME(drot)(&n, dx, &incx, dy, &incy, &c, &s);
}

inline
void drotg(double* restrict da, double* restrict db,
           double* restrict c, double* restrict s)
{
   TRACE_BLAS1("BLAS1: drotg")(da)(db)(c)(s);
   raw::F77NAME(drotg)(da, db, c, s);
}

inline
void dscal(integer n, double da, double* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: dscal")(n)(da)(dx)(incx);
   raw::F77NAME(dscal)(&n, &da, dx, &incx);
}

inline
void dswap(integer n, double* restrict dx, integer incx, double* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: dswap")(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(dswap)(&n, dx, &incx, dy, &incy);
}

inline
integer idamax(integer n, double const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: idmax")(n)(dx)(incx);
   return raw::F77NAME(idamax)(&n, dx, &incx)-1;
}

// complex

inline
double dzasum(integer n, std::complex<double> const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: dzasum")(n)(dx)(incx);
   return raw::F77NAME(dzasum)(&n, dx, &incx);
}

inline
void zaxpy(integer n, std::complex<double> da,
           std::complex<double> const* dx, integer incx,
           std::complex<double>* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: zaxpy")(n)(da)(dx)(incx)(dy)(incy);
   raw::F77NAME(zaxpy)(&n, &da, dx, &incx, dy, &incy);
}

inline
void zcopy(integer n, std::complex<double> const* dx, integer incx,
           std::complex<double>* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: zcopy")(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(zcopy)(&n, dx, &incx, dy, &incy);
}

#if defined(FORTRAN_COMPLEX_RETURN_FIRST_ARG)
inline
std::complex<double>
zdotu(integer n, std::complex<double> const* dx, integer incx,
      std::complex<double> const* dy, integer incy)
{
   std::complex<double> result;
   TRACE_BLAS1("BLAS1: zdotu")(&result)(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(zdotu)(&result, &n, dx, &incx, dy, &incy);
   return result;
}

inline
std::complex<double>
zdotc(integer n, std::complex<double> const* dx, integer incx,
      std::complex<double> const* dy, integer incy)
{
   std::complex<double> result;
   TRACE_BLAS1("BLAS1: zdotc")(&result)(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(zdotc)(&result, &n, dx, &incx, dy, &incy);
   return result;
}
#elif defined(FORTRAN_COMPLEX_RETURN_IN_REGISTER)
inline
std::complex<double>
zdotu(integer n, std::complex<double> const* dx, integer incx,
      std::complex<double> const* dy, integer incy)
{
   TRACE_BLAS1("BLAS1: zdotu")(n)(dx)(incx)(dy)(incy);
   return raw::F77NAME(zdotu)(&n, dx, &incx, dy, &incy);
}

inline
std::complex<double>
zdotc(integer n, std::complex<double> const* dx, integer incx,
      std::complex<double> const* dy, integer incy)
{
   TRACE_BLAS1("BLAS1: zdotc")(&result)(n)(dx)(incx)(dy)(incy);
   return raw::F77NAME(zdotc)(&n, dx, &incx, dy, &incy);
}
#endif

inline
double dznrm2(integer n, std::complex<double> const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: dznrm2")(n)(dx)(incx);
   return raw::F77NAME(dznrm2)(&n, dx, &incx);
}

inline
void zrot(integer n, std::complex<double>* restrict dx, integer incx,
          std::complex<double>* restrict dy, integer incy,
          std::complex<double> c, std::complex<double> s)
{
   TRACE_BLAS1("BLAS1: zrot")(n)(dx)(incx)(dy)(incy)(c)(s);
   raw::F77NAME(zrot)(&n, dx, &incx, dy, &incy, &c, &s);
}

inline
void zrotg(std::complex<double>* restrict da, std::complex<double>* restrict db,
           std::complex<double>* restrict c, std::complex<double>* restrict s)
{
   TRACE_BLAS1("BLAS1: zrotg")(da)(db)(c)(s);
   raw::F77NAME(zrotg)(da, db, c, s);
}

inline
void zscal(integer n, std::complex<double> da, std::complex<double>* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: zscal")(n)(da)(dx)(incx);
   raw::F77NAME(zscal)(&n, &da, dx, &incx);
}

inline
void zdscal(integer n, double da, std::complex<double>* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: zdscal")(n)(da)(dx)(incx);
   raw::F77NAME(zdscal)(&n, &da, dx, &incx);
}

inline
void zswap(integer n, std::complex<double>* restrict dx, integer incx,
           std::complex<double>* restrict dy, integer incy)
{
   TRACE_BLAS1("BLAS1: zswap")(n)(dx)(incx)(dy)(incy);
   raw::F77NAME(zswap)(&n, dx, &incx, dy, &incy);
}


inline
integer izamax(integer n, std::complex<double> const* dx, integer incx)
{
   TRACE_BLAS1("BLAS1: izmax")(n)(dx)(incx);
   return raw::F77NAME(izamax)(&n, dx, &incx)-1;
}

} // namespace BLAS


/* BLAS documentation

**************************************************************************************

sasum, dasum, scasum, dzasum
Sum of the absolute value

FORMAT
  {S,D}ASUM (n,x,incx)
  SCASUM (n,x,incx)
  DZASUM (n,x,incx)

Function Value
  sum: real*4 | real*8 | complex*8 | complex*16
  The sum of the absolute values of the elements of the vector x.
  If n<=0, sum returns the value 0.0.

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

Description
  The SASUM and DASUM functions compute the sum of the absolute values of the
  elements of a real vector x: SUM(i=1...n,|x(i)|) = |x(1)| + |x(2)| + ... +
  |x(n)|

  SCASUM and DZASUM compute the sum of the absolute values of the real and
  imaginary parts of the elements of a complex vector x: SUM(i=1...(n),|a(i)|
  + |b(i)|) = (|a(1)| + |b(1)|) + (|a(2)| + |b(2)|) + ...  + (|a(n)| +
  |b(n)|)

  where x(i) = (a(i),b(i)) and |x(i)| = |a(i)| + |b(i)| = |real| + |ima-
  ginary|

  If incx < 0, the result is identical to using |incx|.  If (incx = 0, the
  computation is a time-consuming way of setting sum = n*x(1).

  Because of the efficient coding of these routines, rounding errors can
  cause the final result to differ from the result computed by a sequential
  evaluation of the sum of the elements of the vector.



Example

  INTEGER*4 N, INCX
  REAL*4 X(20), SUM
  INCX = 1
  N = 20
  SUM = SASUM(N,X,INCX)

  This FORTRAN code shows how to compute the sum of the absolute values of
  the elements of the vector x.


**************************************************************************************


saxpy, daxpy, caxpy,zaxpy
Vector plus the product of a scalar and a vector

FORMAT
  {S,D,C,Z}AXPY (n, alpha, x, incx, y, incy)

Arguments

  n                   integer*4
                      On entry, the number of elements in the vectors x and
                      y.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar multiplier alpha for the elements
                      of the vector x.
                      On exit, lpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|), containing the elements of the vector
                      y.

                      On exit, if n<=0 or alpha = 0, y is unchanged. If n>0,
                      y is overwritten; y(i) is replaced by y(i)+alpha*x(i).

  incy                integer*4
                      On entry, the increment for the array Y.
                      If incy > 0, vector y is stored forward in the array,
                      so that y(i) is stored in location Y(1+(i-1)*incy).
                      If incy < 0, vector y is stored backward in the array,
                      so that y(i) is stored in location Y(1+(n-i)*|incy|).
                      On exit, incy is unchanged.

Description
  The _AXPY functions compute the following scalar-vector product and sum: y
  = alpha*x+y
  where alpha is a scalar, and x and y are vectors.

  If any element of x or the scalar alpha share a memory location with an
  element of y, the results are unpredictable.

  If incx = 0, the computation is a time-consuming way of adding the constant
  alpha*x(1) to all the elements of y.

Example

  INTEGER*4 N, INCX, INCY
  REAL*4 X(20), Y(20), alpha
  INCX = 1
  INCY = 1
  alpha = 2.0
  N = 20
  CALL SAXPY(N,alpha,X,INCX,Y,INCY)

  This FORTRAN code shows how all elements of the real vector x are multi-
  plied by 2.0, added to the elements of the real vector y, and the vector y
  is set equal to the result.


**************************************************************************************


scopy, dcopy, ccopy, zcopy
Copy of a vector

FORMAT
  {S,D,C,Z}COPY (n, x, incx, y, incy)

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).
                      On exit, if n<=0, y is unchanged.  If n > 0, y is
                      overwritten; y(i) is replaced by x(i).

  incy                integer*4
                      On entry, the increment for the array Y.
                      If incy >= 0, vector y is stored forward in the array,
                      so that y(i) is stored in location Y(1+(i-1)*incy).
                      If incy < 0, vector y is stored backward in the array,
                      so that y(i) is stored in location Y(1+(n-i)*|incy|).
                      On exit, incy is unchanged.

Description
  The _COPY subprograms copy the elements of the vector x to the vector y,
  performing the following operation: y(i) = x(i)
  If incx = 0, each y(i) is set to x(1).  Therefore, you can use incx = 0 to
  initialize all elements to a constant.

  If incy = 0, the computation is a time-consuming way of setting y(1) =
  x(n), the last referenced element of the vector x.

  If incy = -incx, the vector x is stored in reverse order in y.  In this
  case, the call format is as follows:

  CALL SCOPY (N,X,INCX,Y,-INCX)


  If any element of x shares a memory location with an element of y, the
  results are unpredictable, except for the following special case.  It is
  possible to move the contents of a vector up or down within itself and not
  cause unpredictable results even though the same memory location is shared
  between input and output. To do this when i > j, call the subroutine with
  incx = incy > 0 as follows:

  CALL SCOPY (N,X(I),INCX,X(J),INCX)


  The call to SCOPY moves elements of the array X x(i),x(i+1*incx),
  ...,x(i+(n-1)*incx) to new elements of the array X x(j),x(j+1*incx), ...,
  x(j+(n-1)*incx). If i < j, specify a negative value for incx and incy in
  the call to the subroutine, as follows. The parts that do not overlap are
  unchanged.

  CALL SCOPY (N,X(I),-INCX,X(J),-INCX)


Examples

  INTEGER*4 N, INCX, INCY
  REAL*4 X(20), Y(20)
  INCX = 1
  INCY = 1
  N = 20
  CALL SCOPY(N,X,INCX,Y,INCY)

  The preceding FORTRAN code copies a vector x to a vector y.

  CALL SCOPY(N,X,-2,X(3),-2))

  The preceding call moves the contents of X(1),X(3),X(5),
   ...  , X(2N-1) to X(3),X(5),
   ...  , X(2N+1) and leaves the vector x unchanged.

  CALL SCOPY(99,X(2),1,X,1))

  The preceding call moves the contents of X(2),X(3),
   ...  , X(100) to X(1),X(2),
   ...  , X(99) and leaves x(100) unchanged.

  CALL SCOPY(N,X,1,Y,-1))

  The preceding call moves the contents of X(1),X(2),X(3),
   ...  , X(N) to Y(N),Y(N-1),
   ...  , Y.


**************************************************************************************


sdot, ddot, dsdot, cdotc, zdotc, cdotu, zdotu
inner product of two vectors

FORMAT
  {S,D}DOT (n, x,incx, y, incy)
  DSDOT (n, x, incx, y, incy)
  {C,Z}DOT{C,U} (n, x, incx, y, incy)

Function Value
  dotpr: real*4 | real*8 | complex*8 | complex*16
  The dot product of the two vectors x and y.

       For real vectors, if n <= 0 , dotpr returns the value 0.0.
       For complex vectors, if n <= 0 , dotpr returns (0.0, 0.0).

Arguments

  n                   integer*4
                      On entry, the number of elements in the vectors x and
                      y.
                      On exit, n is unchanged.

  X                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      y.
                      On exit, y is unchanged.

  incy                integer*4
                      On entry, the increment for the array Y.
                      If incy >= 0, vector y is stored forward in the array,
                      so that y(i) is stored in location Y(1+(i-1)*incy).
                      If incy < 0, vector y is stored backward in the array,
                      so that y(i) is stored in location Y(1+(n-i)*|incy|).
                      On exit, incy is unchanged.




Description
  SDOT, DDOT, and DSDOT compute the dot product of two real vectors.  CDOTC
  and ZDOTC compute the conjugated dot product of two complex vectors.  CDOTU
  and ZDOTU compute the unconjugated dot product of two complex vectors.


  SDOT, DDOT, DSDOT are functions that compute the dot product of two n-
  element real vectors, x and y:
  x dot y = SUM(i=1...n,x(i)y(i)) = x(1)y(1) + x(2)y(2) + ... + x(n)y(n)

  The order of operations is different from the order in a sequential evalua-
  tion of the dot product. The final result can differ from the result of a
  sequential evaluation. The DSDOT functions returns the value in double-
  precision.

  CDOTC and ZDOTC are functions that compute the conjugated dot product of
  two complex vectors, x and y, that is, the complex conjugate of the first
  vector is used to compute the dot product.

  Each element x(j) of the vector x is a complex number and each element y(j)
  of the vector y is a complex number.  The conjugated dot product of two
  complex vectors, x and y, is expressed as follows:
  conjugate(x) dot y = SUM(i=1...n,conjugate(x(i))y(i)) =
  = conjugate(x)(1)y(1) + conjugate(x)(2)y(2) + ... + conjugate(x)(n)y(n)

  For example, x and y each have two complex elements:
  x = (1 + i, 2 - i), y = (3 + i, 3 + 2i)
  The conjugate of vector x is
  conjugate(x) = (1 - i, 2 + i), and the dot product is
  conjugate(x) dot y = (1-i)(3+i) + (2+i)(3+2i) = (4-2i) + (4+7i) = (8+5i))

  CDOTU and ZDOTU compute the unconjugated dot product of two complex vec-
  tors. The unconjugated dot product of two complex vectors, x and y, is
  expressed as follows:
  x dot y = SUM(i=1...n,x(i)y(i)) = x(1)y(1) + x(2)y(2) + ... + x(n)y(n)

  For example, for the same complex vectors x and y:
  x dot y = (1+i)(2+i) + (2-i)(3+2i) = (1+3i) + (8+i) = 9+4i

Example

  INTEGER*4 INCX, INCY
  REAL*4 X(20), Y(20), DOTPR
  INCX = 1
  INCY = 1
  N = 20
  DOTPR = SDOT(N,X,INCX,Y,INCY)

  This FORTRAN code shows how to compute the dot product of two vectors, x
  and y, and return the result in dotpr.

  INTEGER*4 INCX, INCY
  COMPLEX*8 X(20), Y(20), DOTPR
  INCX = 1
  INCY = 1
  N = 20
  DOTPR = CDOTU(N,X,INCX,Y,INCY)

  This FORTRAN code shows how to compute the unconjugated dot product of two
  complex vectors, x and y, and return the result in dotpr.



**************************************************************************************


snrm2, dnrm2, scnrm2, dznrm2
Square root of sum of the squares of the elements of a vector

FORMAT
  {S,D}NRM2 (n, x, incx)
  SCNRM2 (n, x, incx)
  DZNRM2 (n, x, incx)

Function Value
  e_norm: real*4 | real*8
  The Euclidean norm of the vector x, that is, the square root of the conju-
  gated dot product of x with itself.

       If n<=0, e_norm returns the value 0.0.

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

Description
  SNRM2 and DNRM2 compute the Euclidean norm of a real vector; SCNRM2 and
  DZNRM2 compute the Euclidean norm of a complex vector. The Euclidean norm
  is the square root of the conjugated dot product of a vector with itself.

  For real vectors: (SUM(i=1...n,x(i)**(2))**(1/2) = (x(1)**(2) + x(2)**(2)
   + ... + x(n)**(2))**(1/2)
  For complex vectors: (SUM(i=1...n,conjugate(x(i))*x(i))**(1/2) =
  ((conjugate(x)(1) * x(1)) + (conjugate(x)(2) * x(2)) + ... +
  (conjugate(x)(n) * x(n)))**(1/2)

  The order of operations is different from the order in a sequential evalua-
  tion of the Euclidean norm. The final result can differ from the result of
  a sequential evaluation.

  If incx < 0, the result is identical to using |incx|.  If incx = 0, the
  computation is a time-consuming way of setting e_norm =
  (n*x(1)**(2))**(1/2).

Example

  INTEGER*4 INCX, N
  REAL*4 X(20), E_NORM
  INCX = 1
  N = 20
  E_NORM = SNRM2(N,X,INCX)

  This FORTRAN code shows how to compute the Euclidean norm of a real vector.



**************************************************************************************


srot, drot, crot, zrot, csrot, zdrot
Apply givens plane rotation

FORMAT
  {S,D,C,Z}ROT (n, x, incx, y, incy, c, s)
  CSROT (n, x,  incx, y, incy, c, s)
  ZDROT (n, x, incx, y, incy, c, s)

Arguments

  n                   integer*4
                      On entry, the number of elements in the vectors x and
                      y.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, if n<=0 or if c is 1.0 and s is 0.0, x is
                      unchanged.  Otherwise, x is overwritten; X contains the
                      rotated vector x.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).  Y contains the n elements of the
                      vector y.
                      On exit, if n<=0 or if c is 1.0 and s is 0.0, y is
                      unchanged.  Otherwise, y is overwritten; Y contains the
                      rotated vector y.

  incy                integer*4
                      On entry, the increment for the array Y.
                      If incy >= 0, vector y is stored forward in the array,
                      so that y(i) is stored in location Y(1+(i-1)*incy).
                      If incy < 0, vector y is stored backward in the array,
                      so that y(i) is stored in location Y(1+(n-i)*|incy|).
                      On exit, incy is unchanged.

  c                   real*4 | real*8
                      On entry, the first rotation element, that is, the
                      cosine of the angle of rotation.  The argument c is the
                      first rotation element generated by the _ROTG subrou-
                      tines.
                      On exit, c is unchanged.

  s                   real*4 | real*8 | complex*8 | complex*16
                      On entry, the second rotation element, that is, the
                      sine of the angle of rotation.  The argument s is the
                      second rotation element generated by the _ROTG subrou-
                      tines.
                      On exit, s is unchanged.

Description
  SROT and DROT apply a real Givens plane rotation to each element in the
  pair of real vectors, x and y. CSROT and ZDROT apply a real Givens plane
  rotation to elements in the complex vectors, x and y.  CROT and ZROT apply
  a complex Givens plane rotation to each element in the pair of complex vec-
  tors x and y.

  The cosine and sine of the angle of rotation are c and s, respectively, and
  are provided by the BLAS Level 1 _ROTG subroutines.

  The Givens plane rotation for SROT, DROT, CSROT, and ZDROT follows: x(i) =
  c*x(i) + s*y(i) y(i) = -s*x(i) + c*y(i)

  The elements of the rotated vector x are x(i)  = cx(i) + sy(i).
  The elements of the rotated vector y are y(i)  =  -sx(i) + cy(i).

  The Givens plane rotation for CROT and ZROT follows: x(i) = c*x(i) + s*y(i)
  y(i) = -conjugate(s)*x(i) + c*y(i)

  The elements of the rotated vector x are x(i)  = cx(i) + sy(i).
  The elements of the rotated vector y are y(i)  =  -conjugate(s)x(i) +
  cy(i).

  If n<=0 or if c = 1.0 and s = 0.0, x and y are unchanged.  If any element
  of x shares a memory location with an element of y, the results are
  unpredictable.

  These subroutines can be used to introduce zeros selectively into a matrix.

Example

  INTEGER*4 INCX, N
  REAL X(20,20), A, B, C, S
  INCX = 20
  N = 20
  A = X(1,1)
  B = X(2,1)
  CALL SROTG(A,B,C,S)
  CALL SROT(N,X,INCX,X(2,1),INCX,C,S)

  This FORTRAN code shows how to rotate the first two rows of a matrix and
  zero out the element in the first column of the second row.


**************************************************************************************


srotg, drotg, crotg, zrotg
Generate elements for a givens plane rotation

FORMAT
  {S,D,C.Z}ROTG (a, b, c, s)

Arguments

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, the first element of the input vector.
                      On exit, a is overwritten with the rotated element r.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, the second element of the input vector.  On
                      exit, for SROTG and DROTG, b is overwritten with the
                      reconstruction element z.  For CROTG and ZROTG, b is
                      unchanged.

  c                   real*4 | real*8
                      On entry, an unspecified variable.
                      On exit, c is overwritten with the first rotation ele-
                      ment, that is, the cosine of the angle of rotation.

  s                   real*4 | real*8 | complex*8 | complex*16
                      On entry, an unspecified variable.
                      On exit, s is overwritten with the second rotation ele-
                      ment, that is, the sine of the angle of rotation.

Description
  The _ROTG subroutines construct a Givens plane rotation that eliminates the
  second element of a two-element vector and can be used to introduce zeros
  selectively into a matrix.

  Using a and b to represent elements of an input real vector, the SROTG and
  DROTG functions calculate the elements c and s of  an orthogonal matrix
  such that:

   c*a + s*b = r
  -s*a + c*b = 0

  Using a and b to represent elements of an input complex vector, the CROTG
  and ZROTG functions calculate the elements real c and complex s of an
  orthogonal matrix such that:

              c*a + s*b = r
  -conjugate(s)*a + c*b = 0

  A real Givens plane rotation is constructed for values a and b by computing
  values for r, c, s, and z, as follows:

  r=p * (a**(2)+b**(2))**(1/2)

  p = SIGN(a)  if  |a| > |b|
  p = SIGN(b)  if  |a|<=|b|

  c = a/r if r is not equal to 0
  c = 1 if r = 0

  s = b/r if r is not equal to 0
  s = 0 if r = 0

  z = s if |a| > |b|
  z = 1/c if |a|<=|b|, c is not equal to 0, and r is not equal to 0.
  z = 1 if |a|<=|b|, c = 0, and r is not equal to 0.
  z = 0 if r = 0

  SROTG and DROTG can use the reconstruction element z to store the rotation
  elements for future use. The quantities c and s are reconstructed from z as
  follows:

  For |z| = 1, c = 0.0  and  s = 1.0

  For |z| < 1, c = (1-z**(2))**(1/2) and s = z

  For |z| > 1, c = 1/z and s = (1-c**(2))**(1/2)

  A complex Givens plane rotation is constructed for values a and b by com-
  puting values for real c, complex s and complex r, as follows:

  p=(|a|**(2)+|b|**(2))**(1/2)

  q = a/|a|

  r = qp  if  |a| is not equal to 0.
  r = b  if  |a| is equal to 0.

  c = |a|/p if |a| is not equal to 0
  c = 0 if |a| is equal to 0

  s = q*conjugate(b)/p if |a| is not equal to 0
  s = (1.0,0.0) if |a| is equal to 0

  The absolute value used in the above definitions corresponds to the strict
  definition of the absolute value of a complex number.

  The arguments c and s are passed to the _ROT subroutines.

Example

  REAL*4 A, B, C, S
  CALL SROTG(A,B,C,S)

  This FORTRAN code shows how to generate the rotation elements for a vector
  of elements a and b.


**************************************************************************************


sscal, dscal, cscal, zscal, csscal, zdscal
Product of a scalar and a vector

FORMAT
  {S,D,C,Z}SCAL (n, alpha, x, incx)
  CSSCAL (n, alpha, x, incx)
  ZDSCAL (n, alpha, x, incx)

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar value used to multiply the ele-
                      ments of vector x.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, if n<=0 or alpha = 1.0, then x is unchanged.
                      Otherwise, x is overwritten; x(i) is replaced by
                      alpha*x(i).

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx > 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      If incx = 0, only the first element in the array is
                      scaled.
                      On exit, incx is unchanged.

Description
  These routines perform the following operation: x = alpha*x
  SSCAL and DSCAL scale the elements of a real vector by computing the pro-
  duct of the vector and a real scalar alpha*x.  CSCAL and ZSCAL scale the
  elements of a complex vector by computing the product of the vector and a
  complex scalar alpha.  CSSCAL and  ZDCAL scale the elements of a complex
  vector by computing the product of the vector and a real scalar alpha.

  If n<=0 or alpha = 1.0, x is unchanged.

  If incx < 0, the result is identical to using |incx|.

  If alpha = 0.0 or (0.0, 0.0), the computation is a time-consuming way of
  setting all elements of the vector x equal to zero.  Use the BLAS Level 1
  Extensions subroutines _SET to set all the elements of a vector to a
  scalar.

  The _SCAL routines are similar to the BLAS Level 1 Extensions subroutines
  _VCAL routines, but the _VCAL routines use an output vector different from
  the input vector.


Example

  INTEGER*4 INCX, N
  COMPLEX*8 X(20), alpha
  INCX = 1
  alpha = (2.0, 1.0)
  N = 20
  CALL CSCAL(N,alpha,X,INCX)

  This FORTRAN code shows how to scale a complex vector x by the complex
  scalar (2.0, 1.0).


**************************************************************************************


sswap, dswap, cswap, zswap
Exchange the elements of two vectors

FORMAT
  {S,D,C,Z}SWAP (n, x, incx, y, incy)

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, if n<=0, x is unchanged.  If n > 0, x is
                      overwritten; the elements in the array X that are the
                      vector x are overwritten by the vector y.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).
                      On exit, if n<=0, y is unchanged.  If n > 0, y is
                      overwritten; the elements in the array Y that are the
                      vector y are overwritten by the vector x.

  incy                integer*4
                      On entry, the increment for the array Y.
                      If incy >= 0, vector y is stored forward in the array,
                      so that y(i) is stored in location Y(1+(i-1)*incy).
                      If incy < 0, vector y is stored backward in the array,
                      so that y(i) is stored in location Y(1+(n-i)*|incy|).
                      On exit, incy is unchanged.

Description
  These subroutines swap n elements of the vector x with n elements of vector
  y: x<=>y

  If any element of x shares a memory location with an element of y, the
  results are unpredictable.

  If n<=0, x and y are unchanged.


  You can use these subroutines to invert the storage of elements of a vector
  within itself. If incx > 0, each element x(i) is moved from location
  X(1+(i-1)*incx) to location X(1+(n-i)*incx).  The following code fragment
  inverts the storage of elements of a vector within itself:

  NN = N/2
  LHALF = 1+(N-NN)*INCX
  CALL SSWAP(NN,X,INCX,X(LHALF),-INCX)


Example

  INTEGER*4 INCX, INCY, N
  REAL*4 X(20), Y(20)
  INCX = 1
  INCY = 1
  N = 20
  CALL SSWAP(N,X,INCX,Y,INCY)

  The preceding FORTRAN code swaps the contents of vectors x and y.

  INCX = 1
  INCY = -1
  N = 50
  CALL SSWAP(N,X,INCX,X(51),INCY)

  The preceding FORTRAN code inverts the order of storage of the elements of
  x within itself; that is, it moves x(1),...,x(100) to x(100),...,x(1).


**************************************************************************************


isamax, idamax, icamax, izamax
Index of the element of a vector with maximum absolute value

FORMAT
  I{S,D,C,Z}AMAX (n, x, incx)

Function Value
  imax: integer*4
  The index of the element of the vector x that is the largest in absolute
  value of all elements of the vector.  If n<=0, imax returns the value 0.

Arguments

  n                   integer*4
                      On entry, the number of elements in the vector x.
                      On exit, n is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the array X.
                      If incx >= 0, vector x is stored forward in the array,
                      so that x(i) is stored in location X(1+(i-1)*incx).
                      If incx < 0, vector x is stored backward in the array,
                      so that x(i) is stored in location X(1+(n-i)*|incx|).
                      On exit, incx is unchanged.

Description
  These functions determine the first integer i among the elements of the
  vector x such that: |x(i)| = MAX{|x(j)|, j = 1,2, ...,n}
  You can use these functions to obtain the pivots in Gaussian elimination.

  For complex vectors, each element of the vector is a complex number. In
  these subprograms, the absolute value of a complex number is defined as the
  absolute value of the real part plus the absolute value of the imaginary
  part: |x(j)| = |a(j)| + |b(j)| = |real| + |imaginary|

  If incx < 0, the result depends on how the program is processed.  See the
  coding information in this document for a discussion of the possible
  results.  If incx = 0, the computation is a time-consuming way of setting
  imax = 1.

Example

  INTEGER*4 IMAX, N, INCX
  REAL*4 X(40)
  INCX = 2
  N = 20
  IMAX = ISAMAX(N,X,INCX)

  This FORTRAN code shows how to compute the index of a real vector element
  with maximum absolute value.


**************************************************************************************

*/

#endif
