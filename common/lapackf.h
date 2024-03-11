// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/lapackf.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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

//  C++ interface to LAPACK.
//
//  Functions are only getting included as needed - no way near complete!
// The general scheme is that the C++ wrapper functions in the LAPACK namespace do
// some basic mangling/demangling of parameters, such as passing quantities by
// value rather than by pointer where possible, and converting the info parameter
// to a reference.  The wrapper functions forward to the extern "C" functions which
// are in namespace LAPACK::raw, it shouldn't ever be necessary to call these directly.
//
// LAPACK documentation at https://www.netlib.org/lapack/
//
// We could probably use LAPACKE instead of this header,

#if !defined(MPTOOLKIT_COMMON_LAPACK_H)
#define MPTOOLKIT_COMMON_LAPACK_H

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

void dgeqrf(integer m, integer n, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info);

void dgelqf(integer m, integer n, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info);

void dorgqr(integer m, integer n, integer k, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info);

void dorglq(integer m, integer n, integer k, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info);

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

void zgelqf(integer m, integer n, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info);

void zungqr(integer m, integer n, integer k, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info);

void zunglq(integer m, integer n, integer k, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info);

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

void F77NAME(dgeqrf)(integer const* m, integer const* n, double* restrict a, integer const* lda, double* restrict tau, double* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(dgelqf)(integer const* m, integer const* n, double* restrict a, integer const* lda, double* restrict tau,
                     double* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(dorgqr)(integer const* m, integer const* n, integer const* k, double* restrict a, integer const* lda, double* restrict tau, double* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(dorglq)(integer const* m, integer const* n, integer const* k, double* restrict a, integer const* lda, double* restrict tau, double* restrict work, integer const* lwork, integer* restrict info);

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
                     complex* restrict a,
                     integer const* lda, double* restrict s, complex* restrict u,
                     integer const* ldu,
                     complex* restrict vh, integer const* ldvh,
                     complex* restrict work, integer const* lwork,
                     double* restrict rwork, integer* restrict info);

void F77NAME(zhegvx)(integer const* itype, char const* jobz, char const* range, char const* uplo,
                     integer const* n,
                     complex* restrict a, integer const* lda,
                     complex* restrict b, integer const* ldb,
                     double const* vl, double const* vu, integer const* il, integer const* iu,
                     double const* abstol, integer* restrict m, double* restrict w,
                     complex* restrict z,
                     integer const* ldz, complex* restrict work, integer const* lwork,
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

void F77NAME(zgebal)(char const* job, integer const* n, complex* restrict a, integer const* lda, integer* restrict ilo,
                     integer* restrict ihi, double* restrict scale, integer* restrict info);

void F77NAME(zgehrd)(integer const* n, integer const* ilo, integer const* ihi, complex* restrict a, integer const* lda,
                     complex* restrict tau, complex* restrict work, integer* restrict lwork, integer* restrict info);

void F77NAME(zhseqr)(char const* job, char const* compz, integer const* n, integer const* ilo, integer const* ihi,
                     complex* restrict h, integer const* ldh,
                     complex* restrict w, complex* restrict z, integer const* ldz,
                     complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zgeqrf)(integer const* m, integer const* n, complex* restrict a, integer const* lda,
                     complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zgelqf)(integer const* m, integer const* n, complex* restrict a, integer const* lda, complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zungqr)(integer const* m, integer const* n, integer const* k, complex* restrict a, integer const* lda, complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

void F77NAME(zunglq)(integer const* m, integer const* n, integer const* k, complex* restrict a, integer const* lda, complex* restrict tau, complex* restrict work, integer const* lwork, integer* restrict info);

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
void dgeqrf(integer m, integer n, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dgeqrf")(m)(n)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
}

inline
void dgelqf(integer m, integer n, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dgelqf")(m)(n)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(dgelqf)(&m, &n, a, &lda, tau, work, &lwork, &info);
}

inline
void dorgqr(integer m, integer n, integer k, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dorgqr")(m)(n)(k)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

inline
void dorglq(integer m, integer n, integer k, double* restrict a, integer lda, double* restrict tau,
            double* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("dorglq")(m)(n)(k)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(dorglq)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
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
   raw:: F77NAME(zgesvd)(&jobu, &jobvh, &m, &n, reinterpret_cast<complex*>(a), &lda, s, reinterpret_cast<complex*>(u), &ldu,
                         reinterpret_cast<complex*>(vh), &ldvh, reinterpret_cast<complex*>(work), &lwork, rwork, &info);
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
   raw::F77NAME(zhegvx)(&itype, &jobz, &range, &uplo, &n, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(b), &ldb, &vl, &vu, &il, &iu,
                        &abstol, &m, w, reinterpret_cast<complex*>(z), &ldz, reinterpret_cast<complex*>(work), &lwork, rwork, iwork, ifail, &info);
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
            std::complex<double>* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("zgelqf")(m)(n)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zgelqf)(&m, &n, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

void zungqr(integer m, integer n, integer k, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("zungqr")(m)(n)(k)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zungqr)(&m, &n, &k, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

void zunglq(integer m, integer n, integer k, std::complex<double>* restrict a, integer lda, std::complex<double>* restrict tau,
            std::complex<double>* restrict work, integer lwork, integer& restrict info)
{
   TRACE_LAPACK("zunglq")(m)(n)(k)(a)(lda)(tau)(work)(lwork)(info);
   raw::F77NAME(zunglq)(&m, &n, &k, reinterpret_cast<complex*>(a), &lda, reinterpret_cast<complex*>(tau),
                        reinterpret_cast<complex*>(work), &lwork, &info);
}

} // namespce LAPACK

#endif
