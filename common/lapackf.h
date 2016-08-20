// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/lapackf.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
/* -*- C++ -*- $Id$
  lapackf.h

  C++ interface to LAPACK.

  Functions are only getting included as I need them - no way near complete!
*/

#if !defined(LAPACK_H_SDFHJK34893R489UE8ER89FDJOIJ)
#define LAPACK_H_SDFHJK34893R489UE8ER89FDJOIJ

#include "fortran.h"
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include "trace.h"
#include <complex>

#if defined(LAPACK_TRACE_DETAILED)
#define TRACE_LAPACK(Msg) TRACE(Msg)
#else
#define TRACE_LAPACK(Msg) DUMMY_TRACE(Msg)
#endif

namespace LAPACK
{

using namespace Fortran;

// wrapper functions

// real

void dgesv(integer n, integer nrhs, double* a, integer lda, integer* restrict ipiv,
           double* b, integer ldb, integer& restrict info);

void dgels(char trans, integer m, integer n, integer nrhs,
           double* a, integer lda, double* b, integer ldb,
           double* work, integer lwork, integer& restrict info);

void dposv(char uplo, integer n, integer nrhs, double *a, integer lda,
           double *b, integer ldb, integer& restrict info);

void dsyev(char jobz, char uplo, integer n, double* restrict a, integer lda,
           double* restrict w, double* restrict work, integer lwork, integer& restrict info);

void dsyevx(char jobz, char range, char uplo, integer n, double* restrict a, integer lda,
            double vl, double vu, integer il, integer iu, double abstol, int& m,
            double* restrict w, double* restrict z, integer ldz,
            double* restrict work, integer lwork, integer* restrict iwork,
            integer* restrict ifail, integer& restrict info);

void dsyevr(char jobz, char range, char uplo, integer n, double* restrict a, integer lda,
            double vl, double vu, integer il, integer iu, double abstol, int& m,
            double* restrict w, double* restrict z, integer ldz, integer* isuppz,
            double* restrict work,
            integer lwork, integer* restrict iwork, int liwork, integer& restrict info);

void dgesvd(char jobu, char jobvt, integer m, integer n, double* restrict a, integer lda,
            double* restrict s, double* restrict u, integer ldu, double* restrict vt, integer ldvt,
            double* restrict work, integer lwork, integer& restrict info);

void dsygvx(integer itype, char jobz, char range, char uplo, integer n, double* restrict a,
            integer lda, double* restrict b, integer ldb, double vl, double vu, integer il,
            integer iu, double abstol, integer& restrict m, double* restrict w,
            double* restrict z, integer ldz, double* restrict work, integer lwork,
            integer* restrict iwork, integer* restrict ifail, integer& restrict info);

void dpotrf(char uplo, integer n, double* restrict a, integer lda, integer info);

double dlamch(char mach);

// complex

void zposv(char uplo, integer n, integer nrhs, std::complex<double> *a, integer lda,
           std::complex<double> *b, integer ldb, integer& restrict info);

void zheev(char jobz, char uplo, integer n, std::complex<double>* restrict a, integer lda,
           double* restrict w, std::complex<double>* restrict work, integer lwork,
           double* restrict rwork, integer& restrict info);

void zgeev(char jobl, char jobr, integer n, std::complex<double>* restrict a, integer lda,
           std::complex<double>* restrict w, std::complex<double>* restrict vl, integer ldvl,
           std::complex<double>* restrict vr, integer ldvr,
           std::complex<double>* restrict work, integer lwork,
           double* restrict rwork, integer& restrict info);

void zgesvd(char jobu, char jobvh, integer m, integer n, std::complex<double>* restrict a,
            integer lda,
            double* restrict s, std::complex<double>* restrict u, integer ldu,
            std::complex<double>* restrict vh, integer ldvh,
            std::complex<double>* restrict work, integer lwork, double* restrict rwork,
            integer& restrict info);

void zhegvx(integer itype, char jobz, char range, char uplo, integer n,
            std::complex<double>* restrict a, integer lda,
            std::complex<double>* restrict b, integer ldb, double vl, double vu, integer il,
            integer iu,
            double abstol, integer& restrict m, double* restrict w,
            std::complex<double>* restrict z, integer ldz,
            std::complex<double>* restrict work, integer lwork, double* restrict rwork,
            integer* restrict iwork,
            integer* restrict ifail, integer& restrict info);

void zpotrf(char uplo, integer n, std::complex<double>* restrict a, integer lda,
            integer& restrict info);

void zgetrf(integer m, integer n, std::complex<double>* restrict a,
            integer lda, integer* restrict ipiv, integer& restrict info);

void zpotri(char uplo, integer n, std::complex<double>* restrict a, integer lda,
            integer& restrict info);

void zgetri(integer n, std::complex<double>* restrict a, integer lda,
            integer* restrict ipiv, std::complex<double>* restrict work,
            integer lwork, integer& restrict info);

void ztrtri(char uplo, char diag, integer n, std::complex<double>* restrict a,
            integer lda, integer& restrict info);

void zhetrd(char uplo, integer n, std::complex<double>* restrict a, integer lda,
            double* restrict d, double* restrict e, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info);

void zgels(char trans, integer m, integer n, integer nrhs, std::complex<double>* restrict a, integer lda,
           std::complex<double>* restrict b, integer ldb, std::complex<double>* restrict work,
           integer* restrict lwork, integer& info);

void zgeqrf(integer m, integer n, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info);

void zgelqf(integer m, integer n, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info);

// integer ileanv(integer ispec, char const* name, char const* opts, integer n1,
// integer n2 = -1, integer n3 = -1, integer n4 = -1);

void zgebal(char job, integer n, std::complex<double>* restrict a, integer lda, integer& ilo,
                     integer& ihi, double* restrict scale, integer& info);

void zgehrd(integer n, integer ilo, integer ihi, std::complex<double>* restrict a, integer lda,
            std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info);

void zhseqr(char job, char compz, integer n, integer ilo, integer ihi, std::complex<double>* restrict h, integer ldh,
            std::complex<double>* restrict w, std::complex<double>* restrict z, integer ldz,
            std::complex<double>* restrict work, integer lwork, integer& info);

namespace raw
{

extern "C"
{

void F77NAME(dgesv)(integer const* n, integer const* nrhs,
                    double* restrict a, integer const* lda,
                    integer* restrict ipiv,
                    double* restrict b, integer const* ldb, integer* restrict info);

void F77NAME(dgels)(char const* trans, integer const* m, integer const* n,
                    integer const* nrhs, double* restrict a, integer const* lda,
                    double* restrict b, integer const* ldb, double* restrict work,
                    integer const* lwork,
                    integer* restrict info);

void F77NAME(dposv)(char const* uplo, integer const* n, integer const* nrhs,
                    double* restrict a, integer const* lda,
                    double* restrict b, integer const* ldb, integer* restrict info);

void F77NAME(dsyev)(char const* jobz, char const* uplo, integer const* n, double* restrict a,
                    integer const* lda,
                    double* restrict w, double* restrict work, integer const* lwork,
                    integer* restrict info);

void F77NAME(dsyevx)(char const* jobz, char const* range, char const* uplo, integer const* n,
                     double* restrict a, integer const* lda,
                     double const* vl, double const* vu, integer const* il, integer const* iu,
                     double const* abstol, integer* restrict m,
                     double* restrict w, double* restrict z, integer const* ldz,
                     double* restrict work, integer const* lwork, integer* restrict iwork,
                     integer* restrict ifail, integer* restrict info);

void F77NAME(dsyevr)(char const* jobz, char const* range, char const* uplo, integer const* n,
                     double* restrict a,
                     integer const* lda, double const* vl, double const* vu, integer const* il,
                     integer const* iu,
                     double const* abstol, integer* restrict m, double* restrict w,
                     double* restrict z, integer const* ldz, integer* isuppz,
                     double* restrict work, integer const* lwork, integer* restrict iwork,
                     integer const* liwork,
                     integer* info);

void F77NAME(dgesvd)(char const* jobu, char const* jobvt, integer const* m, integer const* n,
                     double* restrict a,
                     integer const* lda, double* restrict s, double* restrict u, integer const* ldu,
                     double* restrict vt, integer const* ldvt, double* restrict work,
                     integer const* lwork,
                     integer* restrict info);

void F77NAME(dsygvx)(integer const* itype, char const* jobz, char const* range, char const* uplo,
                     integer const* n,
                     double* restrict a, integer const* lda, double* restrict b, integer const* ldb,
                     double const* vl, double const* vu, integer const* il, integer const* iu,
                     double const* abstol, integer* restrict m, double* restrict w,
                     double* restrict z, integer const* ldz,
                     double* restrict work, integer const* lwork, integer* restrict iwork,
                     integer* restrict ifail, integer* restrict info);

void F77NAME(dpotrf)(char const* uplo, integer const* n,
                     double* restrict a, integer const* lda, integer* restrict info);

double F77NAME(dlamch)(char const* mach);

void F77NAME(zposv)(char const* uplo, integer const* n, integer const* nrhs, complex* restrict a,
                    integer const* lda,
                    complex* restrict b, integer const* ldb, integer* restrict info);

void F77NAME(zheev)(char const* jobz, char const* uplo, integer const* n, complex* restrict a,
                    integer const* lda,
                    double* restrict w, complex* restrict work, integer const* lwork,
                    double* restrict rwork, integer* restrict info);

void F77NAME(zgeev)(char const* jobl, char const* jobr, integer const* n, complex* restrict a,
                    integer const* lda,
                    complex* restrict w, complex* restrict vl, integer const* ldvl,
                    complex* restrict vr, integer const* ldvr,
                    complex* restrict work, integer const* lwork,
                    double* restrict rwork, integer* restrict info);

void F77NAME(zgesvd)(char const* jobu, char const* jobvh, integer const* m, integer const* n,
                     std::complex<double>* restrict a,
                     integer const* lda, double* restrict s, std::complex<double>* restrict u,
                     integer const* ldu,
                     std::complex<double>* restrict vh, integer const* ldvh,
                     std::complex<double>* restrict work, integer const* lwork,
                     double* restrict rwork, integer* restrict info);

void F77NAME(zhegvx)(integer const* itype, char const* jobz, char const* range, char const* uplo,
                     integer const* n,
                     std::complex<double>* restrict a, integer const* lda,
                     std::complex<double>* restrict b, integer const* ldb,
                     double const* vl, double const* vu, integer const* il, integer const* iu,
                     double const* abstol, integer* restrict m, double* restrict w,
                     std::complex<double>* restrict z,
                     integer const* ldz, std::complex<double>* restrict work, integer const* lwork,
                     double* restrict rwork, integer* restrict iwork, integer* restrict ifail,
                     integer* restrict info);

void F77NAME(zpotrf)(char const* uplo, integer const* n,
                     complex* restrict a, integer const* lda, integer* restrict info);

void F77NAME(zgetrf)(integer const* m, integer const* n, complex* restrict a,
            integer const* lda, integer* restrict ipiv, integer* restrict info);

void F77NAME(zpotri)(char const* uplo, integer const* n,
                     complex* restrict a, integer const* lda, integer* restrict info);

void F77NAME(zgetri)(integer const* n, complex* restrict a, integer const* lda,
            integer* restrict ipiv, complex* restrict work, integer const* lwork,
            integer* restrict info);

void F77NAME(ztrtri)(char const* uplo, char const* diag, integer const* n,
                     complex* restrict a, integer const* lda, integer* restrict info);

void F77NAME(zhetrd)(char const* uplo, integer const* n,
                     complex* restrict a, integer const* lda,
                     double* restrict d, double* restrict e, complex* restrict tau,
                     complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zgels)(char const* trans, integer const* m, integer const* n, integer const* nrhs,
                    complex* restrict a, integer const* lda,
                    complex* restrict b, integer const* ldb, complex* restrict work,
                    integer* restrict lwork, integer* restrict info);

void F77NAME(zgeqrf)(integer const* m, integer const* n, complex* restrict a, integer const* lda,
                     complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zgelqf)(integer const* m, integer const* n, complex* restrict a, integer const* lda,
                     complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zgebal)(char const* job, integer const* n, complex* restrict a, integer const* lda, integer* restrict ilo,
                     integer* restrict ihi, double* restrict scale, integer* restrict info);

void F77NAME(zgehrd)(integer const* n, integer const* ilo, integer const* ihi, complex* restrict a, integer const* lda,
                     complex* restrict tau, complex* restrict work, integer* restrict lwork, integer* restrict info);

void F77NAME(zhseqr)(char const* job, char const* compz, integer const* n, integer const* ilo, integer const* ihi,
                     complex* restrict h, integer const* ldh,
                     complex* restrict w, complex* restrict z, integer const* ldz,
                     complex* restrict work, integer const* lwork, integer* restrict info);


} // extern "C"

} // namespace raw

// implementation of the wrappers

inline
void dgesv(integer n, integer nrhs, double* restrict a, integer lda, integer* restrict ipiv,
           double* restrict b, integer ldb, integer& restrict info)
{
   TRACE_LAPACK("dgesv")(n)(nrhs)(a)(lda)(ipiv)(b)(ldb)(info);
   raw::F77NAME(dgesv)(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}

inline
void dgels(char trans, integer m, integer n, integer nrhs,
           double* a, integer lda, double* b, integer ldb,
           double* work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dgels")(trans)(m)(n)(nrhs)(a)(lda)(b)(ldb)(work)(lwork)(info);
   raw::F77NAME(dgels)(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
}

inline
void dposv(char uplo, integer n, integer nrhs, double* restrict a, integer lda,
           double* restrict b, integer ldb, integer& restrict info)
{
   TRACE_LAPACK("dposv")(uplo)(n)(nrhs)(a)(lda)(b)(ldb)(info);
   raw::F77NAME(dposv)(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
}

inline
void dsyev(char jobz, char uplo, integer n, double* restrict a, integer lda,
           double* restrict w, double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dsyev")(jobz)(uplo)(n)(a)(lda)(w)(work)(lwork)(info);
   raw::F77NAME(dsyev)(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
}

inline
void dsyevx(char jobz, char range, char uplo, integer n, double* restrict a, integer lda,
            double vl, double vu, integer il, integer iu, double abstol, int& m,
            double* restrict w, double* restrict z, integer ldz,
            double* restrict work, integer lwork, integer* restrict iwork,
            integer* restrict ifail, integer& restrict info)
{
   TRACE_LAPACK("dsyevx")(jobz)(range)(uplo)(n)(a)(lda)(vl)(vu)(il)(iu)(abstol)(m)(w)(z)(ldz)(work)(lwork)(iwork)(ifail)(info);
   raw::F77NAME(dsyevx)(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m,
                        w, z, &ldz, work, &lwork, iwork, ifail, &info);
}

inline
void dsyevr(char jobz, char range, char uplo, integer n, double* restrict a, integer lda,
            double vl, double vu, integer il, integer iu, double abstol, int& m,
            double* restrict v, double* restrict z, integer ldz, integer* isuppz,
            double* restrict work,
            integer lwork, integer* restrict iwork, int liwork, integer& restrict info)
{
   TRACE_LAPACK("dsyevr")(jobz)(range)(uplo)(n)(a)(lda)(vl)(vu)(il)(iu)
               (abstol)(m)(v)(z)(ldz)(isuppz)(work)(lwork)(iwork)(liwork)(info);
   raw::F77NAME(dsyevr)(&jobz, &range, &uplo, &n, a, &lda,
                        &vl, &vu, &il, &iu, &abstol, &m,
                        v, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
}

inline
void dgesvd(char jobu, char jobvt, integer m, integer n, double* restrict a, integer lda,
            double* restrict s, double* restrict u, integer ldu, double* restrict vt, integer ldvt,
            double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dgesvd")(jobu)(jobvt)(m)(n)(a)(lda)(s)(u)(ldu)(vt)(ldvt)(work)(lwork)(info);
   raw:: F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
}

inline
void dsygvx(integer itype, char jobz, char range, char uplo, integer n, double* restrict a,
            integer lda,
            double* restrict b, integer ldb, double vl, double vu, integer il, integer iu,
            double abstol, integer& restrict m, double* restrict w, double* restrict z, integer ldz,
            double* restrict work, integer lwork, integer* restrict iwork,
            integer* restrict ifail, integer& restrict info)
{
   TRACE_LAPACK("dsygvx")(itype)(jobz)(range)(uplo)(n)(a)(lda)(b)(ldb)(vl)(vu)(il)(iu)(abstol)
               (m)(w)(z)(ldz)(work)(lwork)(iwork)(ifail)(info);
   raw::F77NAME(dsygvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu,
                        &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
}

inline
void dpotrf(char uplo, integer n, double* restrict a, integer lda, integer info)
{
   TRACE_LAPACK("dpotrf")(uplo)(n)(a)(lda)(info);
   raw::F77NAME(dpotrf)(&uplo, &n, a, &lda, &info);
}

inline
double dlamch(char mach)
{
   return raw::F77NAME(dlamch)(&mach);
}

inline
void zposv(char uplo, integer n, integer nrhs, std::complex<double>* restrict a, integer lda,
           std::complex<double>* restrict b, integer ldb, integer& restrict info)//HERE
{
   TRACE_LAPACK("zposv")(uplo)(n)(nrhs)(a)(lda)(b)(ldb)(info);
   raw::F77NAME(zposv)(&uplo, &n, &nrhs, reinterpret_cast<complex*>(a), &lda,
                       reinterpret_cast<complex*>(b), &ldb, &info);
}

inline
void zheev(char jobz, char uplo, integer n, std::complex<double>* restrict a, integer lda,
           double* restrict w, std::complex<double>* restrict work, integer lwork,
           double* restrict rwork, integer& restrict info)
{
   TRACE_LAPACK("zheev")(jobz)(uplo)(n)(a)(lda)(w)(work)(lwork)(rwork)(info);
   raw::F77NAME(zheev)(&jobz, &uplo, &n, reinterpret_cast<complex*>(a), &lda,
                       w, reinterpret_cast<complex*>(work), &lwork, rwork, &info);
}

inline
void zgeev(char jobl, char jobr, integer n, std::complex<double>* restrict a, integer lda,
           std::complex<double>* restrict w, std::complex<double>* restrict vl, integer ldvl,
           std::complex<double>* restrict vr, integer ldvr,
           std::complex<double>* restrict work, integer lwork,
           double* restrict rwork, integer& restrict info)
{
   TRACE_LAPACK("zgeev")(jobl)(jobr)(n)(a)(lda)(w)(vl)(ldvl)(vr)(ldvr)(work)(lwork)(rwork)(info);
   raw::F77NAME(zgeev)(&jobl, &jobr, &n, reinterpret_cast<complex*>(a), &lda,
                       reinterpret_cast<complex*>(w),
                       reinterpret_cast<complex*>(vl), &ldvl,
                       reinterpret_cast<complex*>(vr), &ldvr,
                       reinterpret_cast<complex*>(work), &lwork, rwork, &info);
}

inline
void zgesvd(char jobu, char jobvh, integer m, integer n, std::complex<double>* restrict a,
            integer lda,
            double* restrict s, std::complex<double>* restrict u, integer ldu,
            std::complex<double>* restrict vh, integer ldvh,
            std::complex<double>* restrict work, integer lwork, double* restrict rwork,
            integer& restrict info)
{
   TRACE_LAPACK("zgesvd")(jobu)(jobvh)(m)(n)(a)(lda)(s)(u)(ldu)(vh)(ldvh)(work)(lwork)(rwork)(info);
   raw:: F77NAME(zgesvd)(&jobu, &jobvh, &m, &n, a, &lda, s, u, &ldu,
                         vh, &ldvh, work, &lwork, rwork, &info);
}


inline
void zhegvx(integer itype, char jobz, char range, char uplo, integer n,
            std::complex<double>* restrict a, integer lda,
            std::complex<double>* restrict b, integer ldb, double vl, double vu, integer il,
            integer iu,
            double abstol, integer& restrict m, double* restrict w,
            std::complex<double>* restrict z, integer ldz,
            std::complex<double>* restrict work, integer lwork, double* restrict rwork,
            integer* restrict iwork,
            integer* restrict ifail, integer& restrict info)
{
   TRACE_LAPACK("zhegvx")(itype)(jobz)(range)(uplo)(n)(a)(lda)(b)(ldb)(vl)(vu)(il)(iu)(abstol)
         (m)(w)(z)(ldz)(work)(lwork)(rwork)(iwork)(ifail)(info);
   raw::F77NAME(zhegvx)(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, &vl, &vu, &il, &iu,
                        &abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, &info);
}

inline
void zpotrf(char uplo, integer n, std::complex<double>* restrict a, integer lda, integer& restrict info)
{
   TRACE_LAPACK("zpotrf")(uplo)(n)(a)(lda)(info);
   raw::F77NAME(zpotrf)(&uplo, &n, reinterpret_cast<complex*>(a), &lda, &info);
}

inline
void zgetrf(integer m, integer n, std::complex<double>* restrict a,
            integer lda, integer* restrict ipiv, integer& restrict info)
{
   TRACE_LAPACK("zgetrf")(m)(n)(a)(lda)(ipiv)(info);
   raw::F77NAME(zgetrf)(&m, &n, reinterpret_cast<complex*>(a), &lda, ipiv, &info);
}

inline
void zpotri(char uplo, integer n, std::complex<double>* restrict a, integer lda,
            integer& restrict info)
{
   TRACE_LAPACK("zpotri")(uplo)(n)(a)(lda)(info);
   raw::F77NAME(zpotri)(&uplo, &n, reinterpret_cast<complex*>(a), &lda, &info);
}

inline
void zgetri(integer n, std::complex<double>* restrict a, integer lda,
            integer* restrict ipiv, std::complex<double>* restrict work,
            integer lwork, integer& restrict info)
{
   TRACE_LAPACK("zpgeri")(n)(a)(lda)(ipiv)(work)(lwork)(info);
   raw::F77NAME(zgetri)(&n, reinterpret_cast<complex*>(a), &lda, ipiv,
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

inline
void ztrtri(char uplo, char diag, integer n, std::complex<double>* restrict a, integer lda,
            integer& restrict info)
{
   TRACE_LAPACK("ztrtri")(uplo)(diag)(n)(a)(lda)(info);
   raw::F77NAME(ztrtri)(&uplo, &diag, &n, reinterpret_cast<complex*>(a), &lda, &info);
}

inline
void zhetrd(char uplo, integer n, std::complex<double>* restrict a, integer lda,
            double* restrict d, double* restrict e, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("zhetrd")(uplo)(n)(a)(lda)(d)(e)(tau)(work)(lwork)(info);
   raw::F77NAME(zhetrd)(&uplo, &n, reinterpret_cast<complex*>(a), &lda,
                        d, e, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

inline
void zgels(char trans, integer m, integer n, integer nrhs, std::complex<double>* restrict a, integer lda,
           std::complex<double>* restrict b, integer ldb, std::complex<double>* restrict work,
           integer* restrict lwork, integer& info)
{
   TRACE_LAPACK("zgels")(trans)(m)(n)(nrhs)(a)(lda)(b)(ldb)(work)(lwork)(info);
   raw::F77NAME(zgels)(&trans, &m, &n, &nrhs, reinterpret_cast<complex*>(a), &lda,
                       reinterpret_cast<complex*>(b), &ldb, reinterpret_cast<complex*>(work), lwork, &info);
}

inline
void zgeqrf(integer m, integer n, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info)
{
   TRACE_LAPACK("zgeqrf")(m)(n)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zgeqrf)(&m, &n, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

inline
void zgelqf(integer m, integer n, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info)
{
   TRACE_LAPACK("zgelqf")(m)(n)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zgelqf)(&m, &n, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

inline
void zgebal(char job, integer n, std::complex<double>* restrict a, integer lda, integer& ilo,
                     integer& ihi, double* restrict scale, integer& info)
{
   TRACE_LAPACK("zgebal")(job)(n)(a)(lda)(ilo)(ihi)(scale)(info);
   raw::F77NAME(zgebal)(&job, &n, reinterpret_cast<complex*>(a), &lda, &ilo,
                        &ihi, scale, &info);
}

inline
void zgehrd(integer n, integer ilo, integer ihi, std::complex<double>* restrict a, integer lda,
            std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& info)
{
   TRACE_LAPACK("zgehrd")(n)(ilo)(ihi)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zgehrd)(&n, &ilo, &ihi, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

inline
void zhseqr(char job, char compz, integer n, integer ilo, integer ihi, std::complex<double>* restrict h, integer ldh,
            std::complex<double>* restrict w, std::complex<double>* restrict z, integer ldz,
            std::complex<double>* restrict work, integer lwork, integer& info)
{
   TRACE_LAPACK("zhseqr")(job)(compz)(n)(ilo)(ihi)(h)(ldh)(w)(z)(ldz)(work)(lwork)(info);
   raw::F77NAME(zhseqr)(&job, &compz, &n, &ilo, &ihi, reinterpret_cast<complex*>(h), &ldh,
                        reinterpret_cast<complex*>(w), reinterpret_cast<complex*>(z), &ldz,
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

} // namespce LAPACK

/*  LAPACK Documentation


      SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DPOSV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N symmetric positive definite matrix and X and B
*  are N-by-NRHS matrices.
*
*  The Cholesky decomposition is used to factor A as
*     A = U**T* U,  if UPLO = 'U', or
*     A = L * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is a lower triangular
*  matrix.  The factored form of A is then used to solve the system of
*  equations A * X = B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i of A is not
*                positive definite, so the factorization could not be
*                completed, and the solution has not been computed.
*
*  =====================================================================


      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )

*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================


     SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
     $                   IWORK, LIWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 20, 2000
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            ISUPPZ( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEVR computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric matrix T.  Eigenvalues and eigenvectors can be
*  selected by specifying either a range of values or a range of
*  indices for the desired eigenvalues.
*
*  Whenever possible, DSYEVR calls DSTEGR to compute the
*  eigenspectrum using Relatively Robust Representations.  DSTEGR
*  computes eigenvalues by the dqds algorithm, while orthogonal
*  eigenvectors are computed from various "good" L D L^T representations
*  (also known as Relatively Robust Representations). Gram-Schmidt
*  orthogonalization is avoided as far as possible. More specifically,
*  the various steps of the algorithm are as follows. For the i-th
*  unreduced block of T,
*     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T
*          is a relatively robust representation,
*     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high
*         relative accuracy by the dqds algorithm,
*     (c) If there is a cluster of close eigenvalues, "choose" sigma_i
*         close to the cluster, and go to step (a),
*     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,
*         compute the corresponding eigenvector by forming a
*         rank-revealing twisted factorization.
*  The desired accuracy of the output can be specified by the input
*  parameter ABSTOL.
*
*  For more details, see "A new O(n^2) algorithm for the symmetric
*  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,
*  Computer Science Division Technical Report No. UCB//CSD-97-971,
*  UC Berkeley, May 1997.
*
*
*  Note 1 : DSYEVR calls DSTEGR when the full spectrum is requested
*  on machines which conform to the ieee-754 floating point standard.
*  DSYEVR calls DSTEBZ and SSTEIN on non-ieee machines and
*  when partial spectrum requests are made.
*
*  Normal execution of DSTEGR may create NaNs and infinities and
*  hence may abort due to a floating point exception in environments
*  which do not handle NaNs and infinities in the ieee standard default
*  manner.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
********** DSTEIN are called
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*          If high relative accuracy is important, set ABSTOL to
*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
*          eigenvalues are computed to high relative accuracy when
*          possible in future releases.  The current code does not
*          make any guarantees about high relative accuracy, but
*          furutre releases will. See J. Barlow and J. Demmel,
*          "Computing Accurate Eigensystems of Scaled Diagonally
*          Dominant Matrices", LAPACK Working Note #7, for a discussion
*          of which matrices define their eigenvalues to high relative
*          accuracy.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, max(1,M))
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ).
********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,26*N).
*          For optimal efficiency, LWORK >= (NB+6)*N,
*          where NB is the max of the blocksize for DSYTRD and DORMTR
*          returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (LIWORK)
*          On exit, if INFO = 0, IWORK(1) returns the optimal LWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  Internal error
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Inderjit Dhillon, IBM Almaden, USA
*     Osni Marques, LBNL/NERSC, USA
*     Ken Stanley, Computer Science Division, University of
*       California at Berkeley, USA
*
* =====================================================================


      SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), S( * ), U( LDU, * ),
     $                   VT( LDVT, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGESVD computes the singular value decomposition (SVD) of a real
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**T, not V.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'A':  all M columns of U are returned in array U:
*          = 'S':  the first min(m,n) columns of U (the left singular
*                  vectors) are returned in the array U;
*          = 'O':  the first min(m,n) columns of U (the left singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
*  JOBVT   (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix
*          V**T:
*          = 'A':  all N rows of V**T are returned in the array VT;
*          = 'S':  the first min(m,n) rows of V**T (the right singular
*                  vectors) are returned in the array VT;
*          = 'O':  the first min(m,n) rows of V**T (the right singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no rows of V**T (no right singular vectors) are
*                  computed.
*
*          JOBVT and JOBU cannot both be 'O'.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          if JOBU = 'O',  A is overwritten with the first min(m,n)
*                          columns of U (the left singular vectors,
*                          stored columnwise);
*          if JOBVT = 'O', A is overwritten with the first min(m,n)
*                          rows of V**T (the right singular vectors,
*                          stored rowwise);
*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
*                          are destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
*          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
*          if JOBU = 'S', U contains the first min(m,n) columns of U
*          (the left singular vectors, stored columnwise);
*          if JOBU = 'N' or 'O', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBU = 'S' or 'A', LDU >= M.
*
*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
*          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
*          V**T;
*          if JOBVT = 'S', VT contains the first min(m,n) rows of
*          V**T (the right singular vectors, stored rowwise);
*          if JOBVT = 'N' or 'O', VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
*          superdiagonal elements of an upper bidiagonal matrix B
*          whose diagonal is in S (not necessarily sorted). B
*          satisfies A = U * B * VT, so it has the same singular values
*          as A, and singular vectors related by U and VT.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 1.
*          LWORK >= MAX(3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if DBDSQR did not converge, INFO specifies how many
*                superdiagonals of an intermediate bidiagonal form B
*                did not converge to zero. See the description of WORK
*                above for details.
*
*  =====================================================================


      SUBROUTINE ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOSV computes the solution to a complex system of linear equations
*     A * X = B,
*  where A is an N-by-N Hermitian positive definite matrix and X and B
*  are N-by-NRHS matrices.
*
*  The Cholesky decomposition is used to factor A as
*     A = U**H* U,  if UPLO = 'U', or
*     A = L * L**H,  if UPLO = 'L',
*  where U is an upper triangular matrix and  L is a lower triangular
*  matrix.  The factored form of A is then used to solve the system of
*  equations A * X = B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**H*U or A = L*L**H.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i of A is not
*                positive definite, so the factorization could not be
*                completed, and the solution has not been computed.
*
*  =====================================================================


      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,2*N-1).
*          For optimal efficiency, LWORK >= (NB+1)*N,
*          where NB is the blocksize for ZHETRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================


      SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGESVD computes the singular value decomposition (SVD) of a complex
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * conjugate-transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**H, not V.
*
*  =====================================================================


      SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGESVD computes the singular value decomposition (SVD) of a complex
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * conjugate-transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**H, not V.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'A':  all M columns of U are returned in array U:
*          = 'S':  the first min(m,n) columns of U (the left singular
*                  vectors) are returned in the array U;
*          = 'O':  the first min(m,n) columns of U (the left singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
*  JOBVT   (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix
*          V**H:
*          = 'A':  all N rows of V**H are returned in the array VT;
*          = 'S':  the first min(m,n) rows of V**H (the right singular
*                  vectors) are returned in the array VT;
*          = 'O':  the first min(m,n) rows of V**H (the right singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no rows of V**H (no right singular vectors) are
*                  computed.
*
*          JOBVT and JOBU cannot both be 'O'.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          if JOBU = 'O',  A is overwritten with the first min(m,n)
*                          columns of U (the left singular vectors,
*                          stored columnwise);
*          if JOBVT = 'O', A is overwritten with the first min(m,n)
*                          rows of V**H (the right singular vectors,
*                          stored rowwise);
*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
*                          are destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (output) COMPLEX*16 array, dimension (LDU,UCOL)
*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
*          If JOBU = 'A', U contains the M-by-M unitary matrix U;
*          if JOBU = 'S', U contains the first min(m,n) columns of U
*          (the left singular vectors, stored columnwise);
*          if JOBU = 'N' or 'O', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBU = 'S' or 'A', LDU >= M.
*
*  VT      (output) COMPLEX*16 array, dimension (LDVT,N)
*          If JOBVT = 'A', VT contains the N-by-N unitary matrix
*          V**H;
*          if JOBVT = 'S', VT contains the first min(m,n) rows of
*          V**H (the right singular vectors, stored rowwise);
*          if JOBVT = 'N' or 'O', VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 1.
*          LWORK >=  2*MIN(M,N)+MAX(M,N).
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (5*min(M,N))
*          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the
*          unconverged superdiagonal elements of an upper bidiagonal
*          matrix B whose diagonal is in S (not necessarily sorted).
*          B satisfies A = U * B * VT, so it has the same singular
*          values as A, and singular vectors related by U and VT.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if ZBDSQR did not converge, INFO specifies how many
*                superdiagonals of an intermediate bidiagonal form B
*                did not converge to zero. See the description of RWORK
*                above for details.
*
*  =====================================================================


      SUBROUTINE SSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     $                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
     $                   LWORK, IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  SSYGVX computes selected eigenvalues, and optionally, eigenvectors
*  of a real generalized symmetric-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
*  and B are assumed to be symmetric and B is also positive definite.
*  Eigenvalues and eigenvectors can be selected by specifying either a
*  range of values or a range of indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A and B are stored;
*          = 'L':  Lower triangle of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrix pencil (A,B).  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) REAL array, dimension (LDA, N)
*          On entry, the symmetric matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*SLAMCH('S').
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          On normal exit, the first M elements contain the selected
*          eigenvalues in ascending order.
*
*  Z       (output) REAL array, dimension (LDZ, max(1,M))
*          If JOBZ = 'N', then Z is not referenced.
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          The eigenvectors are normalized as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*
*          If an eigenvector fails to converge, then that column of Z
*          contains the latest approximation to the eigenvector, and the
*          index of the eigenvector is returned in IFAIL.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) REAL array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,8*N).
*          For optimal efficiency, LWORK >= (NB+3)*N,
*          where NB is the blocksize for SSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  SPOTRF or SSYEVX returned an error code:
*             <= N:  if INFO = i, SSYEVX failed to converge;
*                    i eigenvectors failed to converge.  Their indices
*                    are stored in array IFAIL.
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
* =====================================================================


      SUBROUTINE CHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     $                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
     $                   LWORK, RWORK, IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
      REAL               ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      REAL               RWORK( * ), W( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  CHEGVX computes selected eigenvalues, and optionally, eigenvectors
*  of a complex generalized Hermitian-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
*  B are assumed to be Hermitian and B is also positive definite.
*  Eigenvalues and eigenvectors can be selected by specifying either a
*  range of values or a range of indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit,  the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX array, dimension (LDB, N)
*          On entry, the Hermitian matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**H*U or B = L*L**H.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  VL      (input) REAL
*  VU      (input) REAL
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) REAL
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*SLAMCH('S').
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (output) REAL array, dimension (N)
*          The first M elements contain the selected
*          eigenvalues in ascending order.
*
*  Z       (output) COMPLEX array, dimension (LDZ, max(1,M))
*          If JOBZ = 'N', then Z is not referenced.
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          The eigenvectors are normalized as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*
*          If an eigenvector fails to converge, then that column of Z
*          contains the latest approximation to the eigenvector, and the
*          index of the eigenvector is returned in IFAIL.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) COMPLEX array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,2*N-1).
*          For optimal efficiency, LWORK >= (NB+1)*N,
*          where NB is the blocksize for CHETRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) REAL array, dimension (7*N)
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  CPOTRF or CHEEVX returned an error code:
*             <= N:  if INFO = i, CHEEVX failed to converge;
*                    i eigenvectors failed to converge.  Their indices
*                    are stored in array IFAIL.
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
*  =====================================================================


      SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOTRF computes the Cholesky factorization of a complex Hermitian
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**H * U,  if UPLO = 'U', or
*     A = L  * L**H,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**H*U or A = L*L**H.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  =====================================================================


      SUBROUTINE ZPOTRI( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZPOTRI computes the inverse of a complex Hermitian positive definite
*  matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
*  computed by ZPOTRF.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular factor U or L from the Cholesky
*          factorization A = U**H*U or A = L*L**H, as computed by
*          ZPOTRF.
*          On exit, the upper or lower triangle of the (Hermitian)
*          inverse of A, overwriting the input factor U or L.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the (i,i) element of the factor U or L is
*                zero, and the inverse could not be computed.
*
*  =====================================================================



      SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGELS solves overdetermined or underdetermined complex linear systems
*  involving an M-by-N matrix A, or its conjugate-transpose, using a QR
*  or LQ factorization of A.  It is assumed that A has full rank.
*
*  The following options are provided:
*
*  1. If TRANS = 'N' and m >= n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A*X ||.
*
*  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
*     an underdetermined system A * X = B.
*
*  3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
*     an undetermined system A**H * X = B.
*
*  4. If TRANS = 'C' and m < n:  find the least squares solution of
*     an overdetermined system, i.e., solve the least squares problem
*                  minimize || B - A**H * X ||.
*
*  Several right hand side vectors b and solution vectors x can be
*  handled in a single call; they are stored as the columns of the
*  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
*  matrix X.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          = 'N': the linear system involves A;
*          = 'C': the linear system involves A**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of
*          columns of the matrices B and X. NRHS >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*            if M >= N, A is overwritten by details of its QR
*                       factorization as returned by ZGEQRF;
*            if M <  N, A is overwritten by details of its LQ
*                       factorization as returned by ZGELQF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,NRHS)
*          On entry, the matrix B of right hand side vectors, stored
*          columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
*          if TRANS = 'C'.
*          On exit, if INFO = 0, B is overwritten by the solution
*          vectors, stored columnwise:
*          if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
*          squares solution vectors; the residual sum of squares for the
*          solution in each column is given by the sum of squares of the
*          modulus of elements N+1 to M in that column;
*          if TRANS = 'N' and m < n, rows 1 to N of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'C' and m >= n, rows 1 to M of B contain the
*          minimum norm solution vectors;
*          if TRANS = 'C' and m < n, rows 1 to M of B contain the
*          least squares solution vectors; the residual sum of squares
*          for the solution in each column is given by the sum of
*          squares of the modulus of elements M+1 to N in that column.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= MAX(1,M,N).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          LWORK >= max( 1, MN + max( MN, NRHS ) ).
*          For optimal performance,
*          LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
*          where MN = min(M,N) and NB is the optimum block size.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO =  i, the i-th diagonal element of the
*                triangular factor of A is zero, so that A does not have
*                full rank; the least squares solution could not be
*                computed.
*
*  =====================================================================



      SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQRF computes a QR factorization of a complex M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================


      SUBROUTINE ZGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGELQF computes an LQ factorization of a complex M-by-N matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
*  A(i,i+1:n), and tau in TAU(i).
*
*  =====================================================================


*/

#endif
