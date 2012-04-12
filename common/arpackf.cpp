
#include "arpackf.h"

namespace ARPACK
{

// iparam_t

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

// ipntr_t

void ipntr_t::get_from_raw(int const* raw)
{
   x                    = raw[0];
   y                    = raw[1];
   bx                   = raw[2];
   next_free            = raw[3];
   t                    = raw[4];
   ritz_values          = raw[5];
   ritz_estimates       = raw[6];
   ritz_values_original = raw[7];
   ritz_error_bounds    = raw[8];
   t_eigenvectors       = raw[9];
   np_shifts            = raw[10];
}

void ipntr_t::put_to_raw(int* raw)
{
   raw[0] = x;
   raw[1] = y;
   raw[2] = bx;
   raw[3] = next_free;
   raw[4] = t;
   raw[5] = ritz_values;
   raw[6] = ritz_estimates;
   raw[7] = ritz_values_original;
   raw[8] = ritz_error_bounds;
   raw[9] = t_eigenvectors;
   raw[10] = np_shifts;
}

// implementation of the wrapper functions

void dsaupd(integer* ido, char bmat, integer n, char const* which,
            integer nev, double tol, double* resid,
            integer ncv, double* V, integer ldv,
            iparam_t* iparam, ipntr_t* ipntr, double *workd,
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

void dseupd(bool rvec, char HowMny, logical* select,
            double* d, double* z, integer ldz,
            double sigma, char bmat, integer n,
            char const* which, integer nev, double tol,
            double* resid, integer ncv, double* V,
            integer ldv, iparam_t* iparam, ipntr_t* ipntr,
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

} // namespace ARPACK

