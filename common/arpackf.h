// -*- C++ -*- $Id$
//
// C++ interface to the fortran ARPACK routines.  Only a few wrapper functions are
// implemented so far: dsaupd, dseupd, znaupd, zneupd

#if !defined(ARPACK_H_DSFHFSH389UER89UJU89UER8JJCERJ809)
#define ARPACK_H_DSFHFSH389UER89UJU89UER8JJCERJ809

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if !defined(HAVE_LIBARPACK)
#error "fatal: ARPACK library required to use this header."
#endif

#include "common/fortran.h"

namespace ARPACK
{

using namespace Fortran;

// some structs for the ARPACK parameters

struct iparam_t
{
   integer ishift;
   integer mxiter;
   integer nconv;
   integer mode;
   integer np;
   integer numop;
   integer numopb;
   integer numreo;

   iparam_t() : ishift(1), mxiter(0), nconv(0), mode(1), np(0), numop(0), numopb(0), numreo(0) {}

   void get_from_raw(int const* raw);
   void put_to_raw(int* raw);
};

// pointers into the workd and workl arrays for double precision variants.  
// These are converted between 0-based arrays
// in dn_ipntr_t, and 1-based arrays in the fortran array.
struct dn_ipntr_t
{
   integer x;
   integer y;
   integer bx;
   integer next_free;
   integer t;
   integer ritz_values;
   integer ritz_estimates;
   integer ritz_values_original;
   integer ritz_error_bounds;
   integer t_eigenvectors;
   integer np_shifts;

   dn_ipntr_t() : x(0), y(0), bx(0), next_free(0), t(0), ritz_values(0),
               ritz_estimates(0), ritz_values_original(0), ritz_error_bounds(0),
               t_eigenvectors(0), np_shifts(0) {}

   void get_from_raw(int const* raw);
   void put_to_raw(int* raw);
};

// pointers into the workd and workl arrays for complex variants.  
// These are converted between 0-based arrays
// in zn_ipntr_t, and 1-based arrays in the fortran array.
struct zn_ipntr_t
{
   integer x;
   integer y;
   integer bx;
   integer next_free;
   integer h;
   integer ritz_values;
   integer ritz_estimates;
   integer bounds;
   integer ritz_values_original;
   integer unused_0;
   integer ritz_error_bounds;
   integer h_schur;
   integer h_eigenvectors;
   integer np_shifts;

   zn_ipntr_t() : x(0), y(0), bx(0), next_free(0), h(0), ritz_values(0),
                  ritz_estimates(0), bounds(0),
                  ritz_values_original(0), unused_0(0), ritz_error_bounds(0),
                  h_schur(0), h_eigenvectors(0), np_shifts(0) {}

   void get_from_raw(int const* raw);
   void put_to_raw(int* raw);
};

// wrapper functions to the fortran ARPACK routines

void dsaupd(integer* ido, char bmat, integer n, char const* which,
            integer nev, double tol, double* resid,
            integer ncv, double* V, integer ldv,
            iparam_t* iparam, dn_ipntr_t* ipntr, double *workd,
            double *workl, integer lworkl, integer* info);

void dseupd(bool rvec, char HowMny, logical* select,
            double* d, double* z, integer ldz,
            double sigma, char bmat, integer n,
            char const* which, integer nev, double tol,
            double* resid, integer ncv, double* V,
            integer ldv, iparam_t* iparam, dn_ipntr_t* ipntr,
            double *workd, double *workl,
            integer lworkl, integer* info);

// complex

void znaupd(integer* restrict ido, char bmat, 
            integer n, char which[3],
            integer nev, double tol,
            std::complex<double>* restrict resid,
            integer ncv,
            std::complex<double>* restrict v, integer ldv,
            iparam_t* iparam, zn_ipntr_t ipntr,
            std::complex<double>* restrict workd,
            std::complex<double>* restrict workl,
            integer lworkl,
            double* restrict rwork,
            integer* restrict info);

void zneupd(bool rvec, char howmny,
            logical* restrict select, std::complex<double>* restrict d,
            std::complex<double>* restrict z, integer ldz,
            std::complex<double> sigma,
            std::complex<double>* restrict workev,
            char bmat, 
            integer n, char which[3],
            integer nev, double tol,
            std::complex<double>* restrict resid,
            integer ncv,
            std::complex<double>* restrict v, integer ldv,
            iparam_t* iparam, zn_ipntr_t* ipntr,
            std::complex<double>* restrict workd,
            std::complex<double>* restrict workl,
            integer lworkl,
            double* restrict rwork,
            integer* restrict info);

// the raw functions are in their own namespace, the wrapper functions call these.
namespace raw
{

extern "C"
{

// debug "common" statement.

/*   struct {  */
/*     integer logfil, ndigit, mgetv0; */
/*     integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd; */
/*     integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd; */
/*     integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd; */
/*   } F77NAME(debug); */


// double precision symmetric routines.

void F77NAME(dsaupd)(integer* ido, char const* bmat, integer const* n, char const* which,
                     integer const* nev, double const* tol, double* resid,
                     integer const* ncv, double* V, integer const* ldv,
                     integer* iparam, integer* ipntr, double* workd,
                     double* workl, integer const* lworkl, integer* info);

void F77NAME(dseupd)(logical *rvec, char *HowMny, logical *select,
                     double *d, double *Z, integer *ldz,
                     double *sigma, char *bmat, integer *n,
                     char const* which, integer const* nev, double const* tol,
                     double *resid, integer const* ncv, double *V,
                     integer *ldv, integer *iparam, integer *ipntr,
                     double* workd, double* workl,
                     integer const* lworkl, integer *info);

// double precision nonsymmetric routines.

void F77NAME(dnaupd)(integer *ido, char *bmat, integer *n, char *which,
                     integer *nev, double *tol, double *resid,
                     integer *ncv, double *V, integer *ldv,
                     integer *iparam, integer *ipntr, double *workd,
                     double *workl, integer *lworkl, integer *info);

void F77NAME(dneupd)(logical *rvec, char *HowMny, logical *select,
                     double *dr, double *di, double *Z,
                     integer *ldz, double *sigmar,
                     double *sigmai, double *workev,
                     char *bmat, integer *n, char *which,
                     integer *nev, double *tol, double *resid,
                     integer *ncv, double *V, integer *ldv,
                     integer *iparam, integer *ipntr,
                     double *workd, double *workl,
                     integer *lworkl, integer *info);

// single precision symmetric routines.

void F77NAME(ssaupd)(integer *ido, char *bmat, integer *n, char *which,
                     integer *nev, float *tol, float *resid,
                     integer *ncv, float *V, integer *ldv,
                     integer *iparam, integer *ipntr, float *workd,
                     float *workl, integer *lworkl, integer *info);

void F77NAME(sseupd)(logical *rvec, char *HowMny, logical *select,
                     float *d, float *Z, integer *ldz,
                     float *sigma, char *bmat, integer *n,
                     char *which, integer *nev, float *tol,
                     float *resid, integer *ncv, float *V,
                     integer *ldv, integer *iparam, integer *ipntr,
                     float *workd, float *workl,
                     integer *lworkl, integer *info);

// single precision nonsymmetric routines.

void F77NAME(snaupd)(integer *ido, char *bmat, integer *n, char *which,
                     integer *nev, float *tol, float *resid,
                     integer *ncv, float *V, integer *ldv,
                     integer *iparam, integer *ipntr, float *workd,
                     float *workl, integer *lworkl, integer *info);

void F77NAME(sneupd)(logical *rvec, char *HowMny, logical *select,
                     float *dr, float *di, float *Z,
                     integer *ldz, float *sigmar,
                     float *sigmai, float *workev, char *bmat,
                     integer *n, char *which, integer *nev,
                     float *tol, float *resid, integer *ncv,
                     float *V, integer *ldv, integer *iparam,
                     integer *ipntr, float *workd, float *workl,
                     integer *lworkl, integer *info);

// double precision complex non-hermitian

void F77NAME(znaupd)(integer* restrict ido, char const* bmat, 
                     integer const* n, char const* which,
                     integer const* nev, double const* tol,
                     std::complex<double>* restrict resid,
                     integer const* ncv,
                     std::complex<double>* restrict v, integer const* ldv,
                     integer* restrict iparam, integer* restrict ipntr,
                     std::complex<double>* restrict workd,
                     std::complex<double>* restrict workl,
                     integer const* lworkl,
                     double* restrict rwork,
                     integer* restrict info);

void F77NAME(zneupd)(logical const* rvec, char const* howmny,
                     logical* restrict select, std::complex<double>* restrict d,
                     std::complex<double>* restrict z, integer const* ldz,
                     std::complex<double> const* sigma,
                     std::complex<double>* restrict workev,
                     char const* bmat, 
                     integer const* n, char const* which,
                     integer const* nev, double const* tol,
                     std::complex<double>* restrict resid,
                     integer const* ncv,
                     std::complex<double>* restrict v, integer const* ldv,
                     integer* restrict iparam, integer* restrict ipntr,
                     std::complex<double>* restrict workd,
                     std::complex<double>* restrict workl,
                     integer const* lworkl,
                     double* restrict rwork,
                     integer* restrict info);

} // extern"C"

} // namespace raw

// iparam_t

inline
void iparam_t::get_from_raw(int const* raw)
{
   ishift = raw[0];
   mxiter = raw[2];
   nconv  = raw[4];
   mode   = raw[6];
   np     = raw[7];
   numop  = raw[8];
   numopb = raw[9];
   numreo = raw[10];
}

inline
void iparam_t::put_to_raw(int* raw)
{
   raw[0] = ishift;
   raw[2] = mxiter;
   raw[3] = 1;        // NB must be 1 for current ARPACK
   raw[4] = nconv;
   raw[6] = mode;
   raw[7] = np;
   raw[8] = numop;
   raw[9] = numopb;
   raw[10] = numreo;
}

// dn_ipntr_t

inline
void dn_ipntr_t::get_from_raw(int const* raw)
{
   x                    = raw[0]-1;
   y                    = raw[1]-1;
   bx                   = raw[2]-1;
   next_free            = raw[3]-1;
   t                    = raw[4]-1;
   ritz_values          = raw[5]-1;
   ritz_estimates       = raw[6]-1;
   ritz_values_original = raw[7]-1;
   ritz_error_bounds    = raw[8]-1;
   t_eigenvectors       = raw[9]-1;
   np_shifts            = raw[10]-1;
}

inline
void dn_ipntr_t::put_to_raw(int* raw)
{
   raw[0] = x+1;
   raw[1] = y+1;
   raw[2] = bx+1;
   raw[3] = next_free+1;
   raw[4] = t+1;
   raw[5] = ritz_values+1;
   raw[6] = ritz_estimates+1;
   raw[7] = ritz_values_original+1;
   raw[8] = ritz_error_bounds+1;
   raw[9] = t_eigenvectors+1;
   raw[10] = np_shifts+1;
}

// zn_ipntr_t

inline
void zn_ipntr_t::get_from_raw(int const* raw)
{
   x                    = raw[0]-1;
   y                    = raw[1]-1;
   bx                   = raw[2]-1;
   next_free            = raw[3]-1;
   h                    = raw[4]-1;
   ritz_values          = raw[5]-1;
   ritz_estimates       = raw[6]-1;
   bounds               = raw[7]-1;
   ritz_values_original = raw[8]-1;
   unused_0             = raw[9];
   ritz_error_bounds    = raw[10]-1;
   h_schur              = raw[11]-1;
   h_eigenvectors       = raw[12]-1;
   np_shifts            = raw[13]-1;
}

inline
void zn_ipntr_t::put_to_raw(int* raw)
{
   raw[0] = x+1;
   raw[1] = y+1;
   raw[2] = bx+1;
   raw[3] = next_free+1;
   raw[4] = h+1;
   raw[5] = ritz_values+1;
   raw[6] = ritz_estimates+1;
   raw[7] = bounds+1;
   raw[8] = ritz_values_original+1;
   raw[9] = unused_0;
   raw[10] = ritz_error_bounds+1;
   raw[11] = h_schur+1;
   raw[12] = h_eigenvectors+1;
   raw[13] = np_shifts+1;
}

// implementation of the wrapper functions

inline
void dsaupd(integer* ido, char bmat, integer n, char const* which,
            integer nev, double tol, double* resid,
            integer ncv, double* V, integer ldv,
            iparam_t* iparam, dn_ipntr_t* ipntr, double *workd,
            double *workl, integer lworkl, integer* info)
{
   integer raw_iparam[11];
   integer raw_ipntr[11];

   iparam->put_to_raw(raw_iparam);
   ipntr->put_to_raw(raw_ipntr);

   raw::F77NAME(dsaupd)(ido, &bmat, &n, which,
                        &nev, &tol, resid,
                        &ncv, V, &ldv,
                        raw_iparam, raw_ipntr, workd,
                        workl, &lworkl, info);

   iparam->get_from_raw(raw_iparam);
   ipntr->get_from_raw(raw_ipntr);
}

inline
void dseupd(bool rvec, char HowMny, logical* select,
            double* d, double* z, integer ldz,
            double sigma, char bmat, integer n,
            char const* which, integer nev, double tol,
            double* resid, integer ncv, double* V,
            integer ldv, iparam_t* iparam, dn_ipntr_t* ipntr,
            double *workd, double *workl,
            integer lworkl, integer* info)
{
   logical rvec_raw = rvec;
   integer raw_iparam[11];
   integer raw_ipntr[11];

   iparam->put_to_raw(raw_iparam);
   ipntr->put_to_raw(raw_ipntr);

   raw::F77NAME(dseupd)(&rvec_raw, &HowMny, select,
                        d, z, &ldz,
                        &sigma, &bmat, &n,
                        which, &nev, &tol,
                        resid, &ncv, V,
                        &ldv, raw_iparam, raw_ipntr,
                        workd, workl,
                        &lworkl, info);

   iparam->get_from_raw(raw_iparam);
   ipntr->get_from_raw(raw_ipntr);
}

// complex

inline
void znaupd(integer* restrict ido, char bmat, 
            integer n, char which[3],
            integer nev, double tol,
            std::complex<double>* restrict resid,
            integer ncv,
            std::complex<double>* restrict v, integer ldv,
            iparam_t* iparam, zn_ipntr_t* ipntr,
            std::complex<double>* restrict workd,
            std::complex<double>* restrict workl,
            integer lworkl,
            double* restrict rwork,
            integer* restrict info)
{
   integer raw_iparam[11];
   integer raw_ipntr[11];

   iparam->put_to_raw(raw_iparam);
   ipntr->put_to_raw(raw_ipntr);

   raw::F77NAME(znaupd)(ido, &bmat, &n, which,
                        &nev, &tol, resid,
                        &ncv, v, &ldv,
                        raw_iparam, raw_ipntr, workd,
                        workl, &lworkl, rwork, info);

   iparam->get_from_raw(raw_iparam);
   ipntr->get_from_raw(raw_ipntr);
}

inline
void zneupd(bool rvec, char howmny,
            logical* restrict select, std::complex<double>* restrict d,
            std::complex<double>* restrict z, integer ldz,
            std::complex<double> sigma,
            std::complex<double>* restrict workev,
            char bmat, 
            integer n, char which[3],
            integer nev, double tol,
            std::complex<double>* restrict resid,
            integer ncv,
            std::complex<double>* restrict v, integer ldv,
            iparam_t* iparam, zn_ipntr_t* ipntr,
            std::complex<double>* restrict workd,
            std::complex<double>* restrict workl,
            integer lworkl,
            double* restrict rwork,
            integer* restrict info)
{
   integer raw_iparam[11];
   integer raw_ipntr[11];

   iparam->put_to_raw(raw_iparam);
   ipntr->put_to_raw(raw_ipntr);

   logical irvec = rvec;

   raw::F77NAME(zneupd)(&irvec, &howmny, select, d, z, &ldz, &sigma, workev,
                        &bmat, &n, which,
                        &nev, &tol, resid,
                        &ncv, v, &ldv,
                        raw_iparam, raw_ipntr, workd,
                        workl, &lworkl, rwork, info);

   iparam->get_from_raw(raw_iparam);
   ipntr->get_from_raw(raw_ipntr);
}

} // namespace ARPACK

/*

***************************************************************************
ARPACK documentation

dsaupd
***************************************************************************

c\BeginDoc
c
c\Name: dsaupd
c
c\Description: 
c
c  Reverse communication interface for the Implicitly Restarted Arnoldi 
c  Iteration.  For symmetric problems this reduces to a variant of the Lanczos 
c  method.  This method has been designed to compute approximations to a 
c  few eigenpairs of a linear operator OP that is real and symmetric 
c  with respect to a real positive semi-definite symmetric matrix B, 
c  i.e.
c                   
c       B*OP = (OP')*B.  
c
c  Another way to express this condition is 
c
c       < x,OPy > = < OPx,y >  where < z,w > = z'Bw  .
c  
c  In the standard eigenproblem B is the identity matrix.  
c  ( A' denotes transpose of A)
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  dsaupd is usually called iteratively to solve one of the 
c  following problems:
c
c  Mode 1:  A*x = lambda*x, A symmetric 
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
c           ===> OP = (inv[K - sigma*M])*M  and  B = M. 
c           ===> Shift-and-Invert mode
c
c  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite, 
c           KG symmetric indefinite
c           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
c           ===> Buckling mode
c
c  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
c           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
c           ===> Cayley transformed mode
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call dsaupd 
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first 
c          call to dsaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          dsaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          (If Mode = 2 see remark 5 below)
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3,4 and 5, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c          IDO = 99: done
c          -------------------------------------------------------------
c             
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          Specify which of the Ritz values of OP to compute.
c
c          'LA' - compute the NEV largest (algebraic) eigenvalues.
c          'SA' - compute the NEV smallest (algebraic) eigenvalues.
c          'LM' - compute the NEV largest (in magnitude) eigenvalues.
c          'SM' - compute the NEV smallest (in magnitude) eigenvalues. 
c          'BE' - compute NEV eigenvalues, half from each end of the
c                 spectrum.  When NEV is odd, compute one more from the
c                 high end than from the low end.
c           (see remark 1 below)
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N.
c
c  TOL     Double precision scalar.  (INPUT)
c          Stopping criterion: the relative accuracy of the Ritz value 
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
c          If TOL .LE. 0. is passed a default is set:
c          DEFAULT = DLAMCH('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine DLAMCH).
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT: 
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector. 
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x. (See remark 4 below).
c
c  V       Double precision N by NCV array.  (OUTPUT)
c          The NCV columns of V contain the Lanczos basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to restart
c          the Arnoldi iteration in an implicit fashion.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the current tridiagonal matrix T are returned in
c                      the part of WORKL array corresponding to RITZ.
c                      See remark 6 below.
c          ISHIFT = 1: exact shifts with respect to the reduced 
c                      tridiagonal matrix T.  This is equivalent to 
c                      restarting the iteration with a starting vector 
c                      that is a linear combination of Ritz vectors 
c                      associated with the "wanted" Ritz values.
c          -------------------------------------------------------------
c
c          IPARAM(2) = LEVEC
c          No longer referenced. See remark 2 below.
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed. 
c          On OUTPUT: actual number of Arnoldi update iterations taken. 
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used. 
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3,4,5; See under \Description of dsaupd for the 
c          five modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), dsaupd returns NP, the number
c          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
c          6 below.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.        
c
c  IPNTR   Integer array of length 11.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Lanczos iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in 
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
c          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
c          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
c                    with the Ritz values located in RITZ in WORKL.
c          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
c
c          Note: IPNTR(8:10) is only referenced by dseupd. See Remark 2.
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c          
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD 
c          as temporary workspace during the iteration. Upon termination
c          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
c          subroutine dseupd uses this output.
c          See Data Distribution Note below.  
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least NCV**2 + 8*NCV .
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iterations allowed
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array WORKL is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Informatinal error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -13: NEV and WHICH = 'BE' are incompatable.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization. The user is advised to check that
c                   enough workspace and array storage has been allocated.
c
c
c\Remarks
c  1. The converged Ritz values are always returned in ascending 
c     algebraic order.  The computed Ritz values are approximate
c     eigenvalues of OP.  The selection of WHICH should be made
c     with this in mind when Mode = 3,4,5.  After convergence, 
c     approximate eigenvalues of the original problem may be obtained 
c     with the ARPACK subroutine dseupd. 
c
c  2. If the Ritz vectors corresponding to the converged Ritz values
c     are needed, the user must call dseupd immediately following completion
c     of dsaupd. This is new starting with version 2.1 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL'
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular 
c     linear systems should be solved with L and L' rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L'z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will 
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.   The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c
c  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
c     must do the following. When IDO = 1, Y = OP * X is to be computed.
c     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
c     must overwrite X with A*X. Y is then the solution to the linear set
c     of equations B*Y = A*X.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the 
c     NP = IPARAM(8) shifts in locations: 
c     1   WORKL(IPNTR(11))           
c     2   WORKL(IPNTR(11)+1)         
c                        .           
c                        .           
c                        .      
c     NP  WORKL(IPNTR(11)+NP-1). 
c
c     The eigenvalues of the current tridiagonal matrix are located in 
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
c  DECOMPOSE  D1(N), D2(N,NCV)
c  ALIGN      RESID(I) with D1(I)
c  ALIGN      V(I,J)   with D2(I,J)
c  ALIGN      WORKD(I) with D1(I)     range (1:N)
c  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
c  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
c  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c
c  Cray MPP syntax:
c  ===============
c  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
c  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
c  REPLICATED WORKL(LWORKL)
c  
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
c     Transformations in a k-Step Arnoldi Method". In Preparation.
c
c\Routines called:
c     dsaup2  ARPACK routine that implements the Implicitly Restarted
c             Arnoldi Iteration.
c     dstats  ARPACK routine that initialize timing and other statistics
c             variables.
c     ivout   ARPACK utility routine that prints integers.
c     second  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: saupd.F   SID: 2.7   DATE OF SID: 8/27/96   RELEASE: 2 
c
c\Remarks
c     1. None
c
c\EndLib




***************************************************************************
ARPACK documentation

dseupd
***************************************************************************




c\BeginDoc
c
c\Name: dseupd
c
c\Description: 
c
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) the corresponding approximate eigenvectors,
c
c      (2) an orthonormal (Lanczos) basis for the associated approximate
c          invariant subspace,
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  (Lanczos) basis is always computed.  There is an additional storage cost 
c  of n*nev if both are requested (in this case a separate array Z must be 
c  supplied).
c
c  These quantities are obtained from the Lanczos factorization computed
c  by DSAUPD for the linear operator OP prescribed by the MODE selection
c  (see IPARAM(7) in DSAUPD documentation.)  DSAUPD must be called before
c  this routine is called. These approximate eigenvalues and vectors are 
c  commonly called Ritz values and Ritz vectors respectively.  They are 
c  referred to as such in the comments that follow.   The computed orthonormal 
c  basis for the invariant subspace corresponding to these Ritz values is 
c  referred to as a Lanczos basis.
c
c  See documentation in the header of the subroutine DSAUPD for a definition 
c  of OP as well as other terms and the relation of computed Ritz values 
c  and vectors of OP with respect to the given problem  A*z = lambda*B*z.  
c
c  The approximate eigenvalues of the original problem are returned in
c  ascending algebraic order.  The user may elect to call this routine
c  once for each desired Ritz vector and store it peripherally if desired.
c  There is also the option of computing a selected set of these vectors
c  with a single call.
c
c\Usage:
c  call dseupd 
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
c       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
c
c  RVEC    LOGICAL  (INPUT) 
c          Specifies whether Ritz vectors corresponding to the Ritz value 
c          approximations to the eigenproblem A*z = lambda*B*z are computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors.
c
c  HOWMNY  Character*1  (INPUT) 
c          Specifies how many Ritz vectors are wanted and the form of Z
c          the matrix of Ritz vectors. See remark 1 below.
c          = 'A': compute NEV Ritz vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE.. 
c          If HOWMNY = 'A' , SELECT is used as a workspace for
c          reordering the Ritz values.
c
c  D       Double precision array of dimension NEV.  (OUTPUT)
c          On exit, D contains the Ritz value approximations to the
c          eigenvalues of A*z = lambda*B*z. The values are returned
c          in ascending order. If IPARAM(7) = 3,4,5 then D represents
c          the Ritz values of OP computed by dsaupd transformed to
c          those of the original eigensystem A*z = lambda*B*z. If 
c          IPARAM(7) = 1,2 then the Ritz values of OP are the same 
c          as the those of A*z = lambda*B*z.
c
c  Z       Double precision N by NEV array if HOWMNY = 'A'.  (OUTPUT)
c          On exit, Z contains the B-orthonormal Ritz vectors of the
c          eigensystem A*z = lambda*B*z corresponding to the Ritz
c          value approximations.
c          If  RVEC = .FALSE. then Z is not referenced.
c          NOTE: The array Z may be set equal to first NEV columns of the 
c          Arnoldi/Lanczos basis array V computed by DSAUPD.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ).  In any case,  LDZ .ge. 1.
c
c  SIGMA   Double precision  (INPUT)
c          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
c          IPARAM(7) = 1 or 2.
c
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to DNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, INFO
c
c         must be passed directly to DSEUPD following the last call
c         to DSAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to DSAUPD and the call to DSEUPD.
c
c  Two of these parameters (WORKL, INFO) are also output parameters:
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:4*ncv) contains information obtained in
c          dsaupd.  They are not changed by dseupd.
c          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
c          untransformed Ritz values, the computed error estimates,
c          and the associated eigenvector matrix of H.
c
c          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
c          of the above information computed by dseupd.
c          -------------------------------------------------------------
c          IPNTR(8): pointer to the NCV RITZ values of the original system.
c          IPNTR(9): pointer to the NCV corresponding error bounds.
c          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
c                     of the tridiagonal matrix T. Only referenced by
c                     dseupd if RVEC = .TRUE. See Remarks.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV must be greater than NEV and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from trid. eigenvalue calculation;
c                Information error from LAPACK routine dsteqr.
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4,5.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: NEV and WHICH = 'BE' are incompatible.
c          = -14: DSAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
c          = -16: HOWMNY = 'S' not yet implemented
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
c     1980.
c  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
c     Computer Physics Communications, 53 (1989), pp 169-179.
c  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
c     Implement the Spectral Transformation", Math. Comp., 48 (1987),
c     pp 663-673.
c  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos 
c     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems", 
c     SIAM J. Matr. Anal. Apps.,  January (1993).
c  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
c     for Updating the QR decomposition", ACM TOMS, December 1990,
c     Volume 16 Number 4, pp 369-377.
c
c\Remarks
c  1. The converged Ritz values are always returned in increasing 
c     (algebraic) order.
c
c  2. Currently only HOWMNY = 'A' is implemented. It is included at this
c     stage for the user who wants to incorporate it. 
c
c\Routines called:
c     dsesrt  ARPACK routine that sorts an array X, and applies the
c             corresponding permutation to a matrix A.
c     dsortr  dsortr  ARPACK sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     dvout   ARPACK utility routine that prints vectors.
c     dgeqr2  LAPACK routine that computes the QR factorization of
c             a matrix.
c     dlacpy  LAPACK matrix copy routine.
c     dlamch  LAPACK routine that determines machine constants.
c     dorm2r  LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
c             of a tridiagonal matrix.
c     dger    Level 2 BLAS rank one update to a matrix.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
c     dscal   Level 1 BLAS that scales a vector.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.

c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational & 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c 
c\Revision history:
c     12/15/93: Version ' 2.1'
c
c\SCCS Information: @(#) 
c FILE: seupd.F   SID: 2.8   DATE OF SID: 5/19/98   RELEASE: 2
c
c\EndLib
c




***************************************************************************
ARPACK documentation

znaupd
***************************************************************************


c\BeginDoc
c
c\Name: znaupd
c
c\Description:
c  Reverse communication interface for the Implicitly Restarted Arnoldi
c  iteration. This is intended to be used to find a few eigenpairs of a
c  complex linear operator OP with respect to a semi-inner product defined
c  by a hermitian positive semi-definite real matrix B. B may be the identity
c  matrix.  NOTE: if both OP and B are real, then dsaupd or dnaupd should
c  be used.
c
c
c  The computed approximate eigenvalues are called Ritz values and
c  the corresponding approximate eigenvectors are called Ritz vectors.
c
c  znaupd is usually called iteratively to solve one of the
c  following problems:
c
c  Mode 1:  A*x = lambda*x.
c           ===> OP = A  and  B = I.
c
c  Mode 2:  A*x = lambda*M*x, M hermitian positive definite
c           ===> OP = inv[M]*A  and  B = M.
c           ===> (If M can be factored see remark 3 below)
c
c  Mode 3:  A*x = lambda*M*x, M hermitian semi-definite
c           ===> OP =  inv[A - sigma*M]*M   and  B = M.
c           ===> shift-and-invert mode
c           If OP*x = amu*x, then lambda = sigma + 1/amu.
c
c
c  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
c        should be accomplished either by a direct method
c        using a sparse matrix factorization and solving
c
c           [A - sigma*M]*w = v  or M*w = v,
c
c        or through an iterative method for solving these
c        systems.  If an iterative method is used, the
c        convergence test must be more stringent than
c        the accuracy requirements for the eigenvalue
c        approximations.
c
c\Usage:
c  call znaupd
c     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
c       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.  IDO must be zero on the first
c          call to znaupd.  IDO will be set internally to
c          indicate the type of operation to be performed.  Control is
c          then given back to the calling routine which has the
c          responsibility to carry out the requested operation and call
c          znaupd with the result.  The operand is given in
c          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    This is for the initialization phase to force the
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c                    In mode 3, the vector B * X is already
c                    available in WORKD(ipntr(3)).  It does not
c                    need to be recomputed in forming OP * X.
c          IDO =  2: compute  Y = M * X  where
c                    IPNTR(1) is the pointer into WORKD for X,
c                    IPNTR(2) is the pointer into WORKD for Y.
c          IDO =  3: compute and return the shifts in the first
c                    NP locations of WORKL.
c          IDO = 99: done
c          -------------------------------------------------------------
c          After the initialization phase, when the routine is used in
c          the "shift-and-invert" mode, the vector M * X is already
c          available and does not need to be recomputed in forming OP*X.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of the matrix B that defines the
c          semi-inner product for the operator OP.
c          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
c          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  WHICH   Character*2.  (INPUT)
c          'LM' -> want the NEV eigenvalues of largest magnitude.
c          'SM' -> want the NEV eigenvalues of smallest magnitude.
c          'LR' -> want the NEV eigenvalues of largest real part.
c          'SR' -> want the NEV eigenvalues of smallest real part.
c          'LI' -> want the NEV eigenvalues of largest imaginary part.
c          'SI' -> want the NEV eigenvalues of smallest imaginary part.
c
c  NEV     Integer.  (INPUT)
c          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
c
c  TOL     Double precision  scalar.  (INPUT)
c          Stopping criteria: the relative accuracy of the Ritz value
c          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
c          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
c          DEFAULT = dlamch('EPS')  (machine precision as computed
c                    by the LAPACK auxiliary subroutine dlamch).
c
c  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
c          On INPUT:
c          If INFO .EQ. 0, a random initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          On OUTPUT:
c          RESID contains the final residual vector.
c
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V. NCV must satisfy the two
c          inequalities 1 <= NCV-NEV and NCV <= N.
c          This will indicate how many Arnoldi vectors are generated
c          at each iteration.  After the startup phase in which NEV
c          Arnoldi vectors are generated, the algorithm generates
c          approximately NCV-NEV Arnoldi vectors at each subsequent update
c          iteration. Most of the cost in generating each Arnoldi vector is
c          in the matrix-vector operation OP*x. (See remark 4 below.)
c
c  V       Complex*16 array N by NCV.  (OUTPUT)
c          Contains the final set of Arnoldi basis vectors.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling program.
c
c  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
c          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
c          The shifts selected at each iteration are used to filter out
c          the components of the unwanted eigenvector.
c          -------------------------------------------------------------
c          ISHIFT = 0: the shifts are to be provided by the user via
c                      reverse communication.  The NCV eigenvalues of
c                      the Hessenberg matrix H are returned in the part
c                      of WORKL array corresponding to RITZ.
c          ISHIFT = 1: exact shifts with respect to the current
c                      Hessenberg matrix H.  This is equivalent to
c                      restarting the iteration from the beginning
c                      after updating the starting vector with a linear
c                      combination of Ritz vectors associated with the
c                      "wanted" eigenvalues.
c          ISHIFT = 2: other choice of internal shift to be defined.
c          -------------------------------------------------------------
c
c          IPARAM(2) = No longer referenced
c
c          IPARAM(3) = MXITER
c          On INPUT:  maximum number of Arnoldi update iterations allowed.
c          On OUTPUT: actual number of Arnoldi update iterations taken.
c
c          IPARAM(4) = NB: blocksize to be used in the recurrence.
c          The code currently works only for NB = 1.
c
c          IPARAM(5) = NCONV: number of "converged" Ritz values.
c          This represents the number of Ritz values that satisfy
c          the convergence criterion.
c
c          IPARAM(6) = IUPD
c          No longer referenced. Implicit restarting is ALWAYS used.
c
c          IPARAM(7) = MODE
c          On INPUT determines what type of eigenproblem is being solved.
c          Must be 1,2,3; See under \Description of znaupd for the
c          four modes available.
c
c          IPARAM(8) = NP
c          When ido = 3 and the user provides shifts through reverse
c          communication (IPARAM(1)=0), _naupd returns NP, the number
c          of shifts the user is to provide. 0 < NP < NCV-NEV.
c
c          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
c          OUTPUT: NUMOP  = total number of OP*x operations,
c                  NUMOPB = total number of B*x operations if BMAT='G',
c                  NUMREO = total number of steps of re-orthogonalization.
c
c  IPNTR   Integer array of length 14.  (OUTPUT)
c          Pointer to mark the starting locations in the WORKD and WORKL
c          arrays for matrices/vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X in WORKD.
c          IPNTR(2): pointer to the current result vector Y in WORKD.
c          IPNTR(3): pointer to the vector B * X in WORKD when used in
c                    the shift-and-invert mode.
c          IPNTR(4): pointer to the next available location in WORKL
c                    that is untouched by the program.
c          IPNTR(5): pointer to the NCV by NCV upper Hessenberg
c                    matrix H in WORKL.
c          IPNTR(6): pointer to the  ritz value array  RITZ
c          IPNTR(7): pointer to the (projected) ritz vector array Q
c          IPNTR(8): pointer to the error BOUNDS array in WORKL.
c          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
c
c          Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.
c
c          IPNTR(9): pointer to the NCV RITZ values of the
c                    original system.
c          IPNTR(10): Not Used
c          IPNTR(11): pointer to the NCV corresponding error bounds.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c
c          -------------------------------------------------------------
c
c  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The user should not use WORKD
c          as temporary workspace during the iteration !!!!!!!!!!
c          See Data Distribution Note below.
c
c  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.  See Data Distribution Note below.
c
c  LWORKL  Integer.  (INPUT)
c          LWORKL must be at least 3*NCV**2 + 5*NCV.
c
c  RWORK   Double precision  work array of length NCV (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c
c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the
c                Implicitly restarted Arnoldi iteration. One possibility
c                is to increase the size of NCV relative to NEV.
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   User input error highly likely.  Please
c                   check actual array dimensions and layout.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.
c
c\Remarks
c  1. The computed Ritz values are approximate eigenvalues of OP. The
c     selection of WHICH should be made with this in mind when using
c     Mode = 3.  When operating in Mode = 3 setting WHICH = 'LM' will
c     compute the NEV eigenvalues of the original problem that are
c     closest to the shift SIGMA . After convergence, approximate eigenvalues
c     of the original problem may be obtained with the ARPACK subroutine zneupd.
c
c  2. If a basis for the invariant subspace corresponding to the converged Ritz
c     values is needed, the user must call zneupd immediately following
c     completion of znaupd. This is new starting with release 2 of ARPACK.
c
c  3. If M can be factored into a Cholesky factorization M = LL`
c     then Mode = 2 should not be selected.  Instead one should use
c     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
c     linear systems should be solved with L and L` rather
c     than computing inverses.  After convergence, an approximate
c     eigenvector z of the original problem is recovered by solving
c     L`z = x  where x is a Ritz vector of OP.
c
c  4. At present there is no a-priori analysis to guide the selection
c     of NCV relative to NEV.  The only formal requirement is that NCV > NEV + 1.
c     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
c     the same type are to be solved, one should experiment with increasing
c     NCV while keeping NEV fixed for a given test problem.  This will
c     usually decrease the required number of OP*x operations but it
c     also increases the work and storage required to maintain the orthogonal
c     basis vectors.  The optimal "cross-over" with respect to CPU time
c     is problem dependent and must be determined empirically.
c     See Chapter 8 of Reference 2 for further information.
c
c  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) complex shifts in locations
c     WORKL(IPNTR(14)), WORKL(IPNTR(14)+1), ... , WORKL(IPNTR(14)+NP).
c     Eigenvalues of the current upper Hessenberg matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are ordered
c     according to the order defined by WHICH.  The associated Ritz estimates
c     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... ,
c     WORKL(IPNTR(8)+NCV-1).
c
c-----------------------------------------------------------------------
c
c\Data Distribution Note:
c
c  Fortran-D syntax:
c  ================
c  Complex*16 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
c  decompose  d1(n), d2(n,ncv)
c  align      resid(i) with d1(i)
c  align      v(i,j)   with d2(i,j)
c  align      workd(i) with d1(i)     range (1:n)
c  align      workd(i) with d1(i-n)   range (n+1:2*n)
c  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
c  distribute d1(block), d2(block,:)
c  replicated workl(lworkl)
c
c  Cray MPP syntax:
c  ===============
c  Complex*16 resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
c  shared     resid(block), v(block,:), workd(block,:)
c  replicated workl(lworkl)
c
c  CM2/CM5 syntax:
c  ==============
c




***************************************************************************
ARPACK documentation

zneupd
***************************************************************************


c\BeginDoc
c
c\Name: zneupd
c
c\Description:
c  This subroutine returns the converged approximations to eigenvalues
c  of A*z = lambda*B*z and (optionally):
c
c      (1) The corresponding approximate eigenvectors;
c
c      (2) An orthonormal basis for the associated approximate
c          invariant subspace;
c
c      (3) Both.
c
c  There is negligible additional cost to obtain eigenvectors.  An orthonormal
c  basis is always computed.  There is an additional storage cost of n*nev
c  if both are requested (in this case a separate array Z must be supplied).
c
c  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
c  are derived from approximate eigenvalues and eigenvectors of
c  of the linear operator OP prescribed by the MODE selection in the
c  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
c  These approximate eigenvalues and vectors are commonly called Ritz
c  values and Ritz vectors respectively.  They are referred to as such
c  in the comments that follow.   The computed orthonormal basis for the
c  invariant subspace corresponding to these Ritz values is referred to as a
c  Schur basis.
c
c  The definition of OP as well as other terms and the relation of computed
c  Ritz values and vectors of OP with respect to the given problem
c  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief
c  description, see definitions of IPARAM(7), MODE and WHICH in the
c  documentation of ZNAUPD.
c
c\Usage:
c  call zneupd
c     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT,
c       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD,
c       WORKL, LWORKL, RWORK, INFO )
c
c\Arguments:
c  RVEC    LOGICAL  (INPUT)
c          Specifies whether a basis for the invariant subspace corresponding
c          to the converged Ritz value approximations for the eigenproblem
c          A*z = lambda*B*z is computed.
c
c             RVEC = .FALSE.     Compute Ritz values only.
c
c             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
c                                See Remarks below.
c
c  HOWMNY  Character*1  (INPUT)
c          Specifies the form of the basis for the invariant subspace
c          corresponding to the converged Ritz values that is to be computed.
c
c          = 'A': Compute NEV Ritz vectors;
c          = 'P': Compute NEV Schur vectors;
c          = 'S': compute some of the Ritz vectors, specified
c                 by the logical array SELECT.
c
c  SELECT  Logical array of dimension NCV.  (INPUT)
c          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
c          computed. To select the  Ritz vector corresponding to a
c          Ritz value D(j), SELECT(j) must be set to .TRUE..
c          If HOWMNY = 'A' or 'P', SELECT need not be initialized
c          but it is used as internal workspace.
c
c  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
c          On exit, D contains the  Ritz  approximations
c          to the eigenvalues lambda for A*z = lambda*B*z.
c
c  Z       Complex*16 N by NEV array    (OUTPUT)
c          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
c          Z represents approximate eigenvectors (Ritz vectors) corresponding
c          to the NCONV=IPARAM(5) Ritz values for eigensystem
c          A*z = lambda*B*z.
c
c          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.
c
c          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
c          the array Z may be set equal to first NEV+1 columns of the Arnoldi
c          basis array V computed by ZNAUPD.  In this case the Arnoldi basis
c          will be destroyed and overwritten with the eigenvector basis.
c
c  LDZ     Integer.  (INPUT)
c          The leading dimension of the array Z.  If Ritz vectors are
c          desired, then  LDZ .ge.  max( 1, N ) is required.
c          In any case,  LDZ .ge. 1 is required.
c
c  SIGMA   Complex*16  (INPUT)
c          If IPARAM(7) = 3 then SIGMA represents the shift.
c          Not referenced if IPARAM(7) = 1 or 2.
c
c  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)
c
c  **** The remaining arguments MUST be the same as for the   ****
c  **** call to ZNAUPD that was just completed.               ****
c
c  NOTE: The remaining arguments
c
c           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
c           WORKD, WORKL, LWORKL, RWORK, INFO
c
c         must be passed directly to ZNEUPD following the last call
c         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
c         the the last call to ZNAUPD and the call to ZNEUPD.
c
c  Three of these parameters (V, WORKL and INFO) are also output parameters:
c
c  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)
c
c          Upon INPUT: the NCV columns of V contain the Arnoldi basis
c                      vectors for OP as constructed by ZNAUPD .
c
c          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
c                       contain approximate Schur vectors that span the
c                       desired invariant subspace.
c
c          NOTE: If the array Z has been set equal to first NEV+1 columns
c          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
c          Arnoldi basis held by V has been overwritten by the desired
c          Ritz vectors.  If a separate array Z has been passed then
c          the first NCONV=IPARAM(5) columns of V will contain approximate
c          Schur vectors that span the desired invariant subspace.
c
c  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
c          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
c          znaupd.  They are not changed by zneupd.
c          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
c          untransformed Ritz values, the untransformed error estimates of
c          the Ritz values, the upper triangular matrix for H, and the
c          associated matrix representation of the invariant subspace for H.
c
c          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
c          of the above information computed by zneupd.
c          -------------------------------------------------------------
c          IPNTR(9):  pointer to the NCV RITZ values of the
c                     original system.
c          IPNTR(10): Not used
c          IPNTR(11): pointer to the NCV corresponding error estimates.
c          IPNTR(12): pointer to the NCV by NCV upper triangular
c                     Schur matrix for H.
c          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
c                     of the upper Hessenberg matrix H. Only referenced by
c                     zneupd if RVEC = .TRUE. See Remark 2 below.
c          -------------------------------------------------------------
c
c  INFO    Integer.  (OUTPUT)
c          Error flag on output.
c          =  0: Normal exit.
c
c          =  1: The Schur form computed by LAPACK routine csheqr
c                could not be reordered by LAPACK routine ztrsen.
c                Re-enter subroutine zneupd with IPARAM(5)=NCV and
c                increase the size of the array D to have
c                dimension at least dimension NCV and allocate at least NCV
c                columns for Z. NOTE: Not necessary if Z and V share
c                the same space. Please notify the authors if this error
c                occurs.
c
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 1 and less than or equal to N.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work WORKL array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation.
c                This should never happened.
c          = -9: Error return from calculation of eigenvectors.
c                Informational error from LAPACK routine ztrevc.
c          = -10: IPARAM(7) must be 1,2,3
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
c          = -12: HOWMNY = 'S' not yet implemented
c          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
c          = -14: ZNAUPD did not find any eigenvalues to sufficient
c                 accuracy.
c          = -15: ZNEUPD got a different count of the number of converged
c                 Ritz values than ZNAUPD got.  This indicates the user
c                 probably made an error in passing data from ZNAUPD to
c                 ZNEUPD or that the data was modified before entering
c                 ZNEUPD
c
c\BeginLib
c
c\References:
c  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
c     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
c     pp 357-385.
c  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
c     Restarted Arnoldi Iteration", Rice University Technical Report
c     TR95-13, Department of Computational and Applied Mathematics.
c  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
c     "How to Implement the Spectral Transformation", Math Comp.,
c     Vol. 48, No. 178, April, 1987 pp. 664-673.
c
c\Routines called:
c     ivout   ARPACK utility routine that prints integers.
c     zmout   ARPACK utility routine that prints matrices
c     zvout   ARPACK utility routine that prints vectors.
c     zgeqr2  LAPACK routine that computes the QR factorization of
c             a matrix.
c     zlacpy  LAPACK matrix copy routine.
c     zlahqr  LAPACK routine that computes the Schur form of a
c             upper Hessenberg matrix.
c     zlaset  LAPACK matrix initialization routine.
c     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
c             in upper triangular form.
c     ztrsen  LAPACK routine that re-orders the Schur form.
c     zunm2r  LAPACK routine that applies an orthogonal matrix in
c             factored form.
c     dlamch  LAPACK routine that determines machine constants.
c     ztrmm   Level 3 BLAS matrix times an upper triangular matrix.
c     zgeru   Level 2 BLAS rank one update to a matrix.
c     zcopy   Level 1 BLAS that copies one vector to another .
c     zscal   Level 1 BLAS that scales a vector.
c     zdscal  Level 1 BLAS that scales a complex vector by a real number.
c     dznrm2  Level 1 BLAS that computes the norm of a complex vector.
c
c\Remarks
c
c  1. Currently only HOWMNY = 'A' and 'P' are implemented.
c
c  2. Schur vectors are an orthogonal representation for the basis of
c     Ritz vectors. Thus, their numerical properties are often superior.
c     If RVEC = .true. then the relationship
c             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
c       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
c     are approximately satisfied.
c     Here T is the leading submatrix of order IPARAM(5) of the
c     upper triangular matrix stored workl(ipntr(12)).
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Chao Yang                    Houston, Texas
c     Dept. of Computational &
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: neupd.F   SID: 2.8   DATE OF SID: 07/21/02   RELEASE: 2
c
c\EndLib



*/

#endif
