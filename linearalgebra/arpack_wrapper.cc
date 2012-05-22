// -*- C++ -*- $Id$ 

#include "common/arpackf.h"
#include "linearalgebra/vectormemproxy.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"

namespace LinearAlgebra
{

template <typename MultFunc>
Vector<std::complex<double> > 
DiagonalizeARPACK(MultFunc Mult, int n, int NumEigen, double tol,
                  std::vector<std::complex<double> >* OutputVectors,
                  int ncv, bool Sort, int Verbose)
{
   if (Verbose >= 1)
   {
      std::cerr << "Total dimension = " << n << '\n';
   }

   LinearAlgebra::Vector<std::complex<double> > Result;

   if (NumEigen >= n)
   {
      // This is a case that ARPACK doesn't handle - it is apparently not
      // capable of producting more than n-1 eigenvalues.  Instead we handle this
      // case by converting to a dense matrix
      if (Verbose >= 1)
      {
         std::cerr << "Constructing matrix for direct diagonalization\n";
      }
      LinearAlgebra::Vector<std::complex<double> > e(n, 0.0), Out(n);
      LinearAlgebra::Matrix<std::complex<double> > Mat(n, n);
      for (int k = 0; k < n; ++k)
      {
         e[k] = 1.0;
         Mult(data(e), data(Out));
         Mat(LinearAlgebra::all, k) = Out;
         e[k] = 0.0;
      }
      LinearAlgebra::Matrix<std::complex<double> > LV, RV;
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
   int const nev = std::min(NumEigen, n-1); // number of eigenvalues to be computed
   std::vector<std::complex<double> > resid(n);  // residual
   ncv = std::min(std::max(ncv, 2*nev + 10), n);            // length of the arnoldi sequence
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
 
   if (Verbose >= 1)
   {
      std::cerr << "Starting ARPACK";
   }

   int NumMultiplies = 0;

   //   TRACE(tol);

   ARPACK::znaupd(&ido, bmat, n, which, nev, tol, &resid[0], ncv,
                  &v[0], ldv, &iparam, &ipntr, &workd[0],
                  &workl[0], lworkl, &rwork[0], &info);
   CHECK(info >= 0)(info)(n);
   
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

   Result = LinearAlgebra::Vector<std::complex<double> >(nev);
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

   return Result;
}

} // namespace LinearAlgebra
