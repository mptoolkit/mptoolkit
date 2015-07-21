// -*- C++ -*- $Id$
//
//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no onvergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

#if !defined(GMRES_H_DHFJDKH48Y78932Y78YHFO7H7O8W)
#define GMRES_H_DHFJDKH48Y78932Y78YHFO7H7O8W

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;
using LinearAlgebra::norm_2_sq;
using LinearAlgebra::norm_2;

double const DGKS_Threshold = 1.0 / std::sqrt(2.0); // 1.0; // between 0 and 1.

template <typename Matrix, typename Vector1, typename Vector2>
void 
Update(Vector1& x, int k, Matrix const& h, Vector2 const& s, Vector1 v[])
{
  Vector2 y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h(j,i) * y[i];
  }

  for (int j = 0; j <= k; j++)
    x += v[j] * y[j];
}

template<typename Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
   if (dy == 0.0) 
   {
      cs = 1.0;
      sn = 0.0;
   } 
   else if (norm_2_sq(dy) > norm_2_sq(dx)) 
   {
      Real temp = dx / dy;
      sn = 1.0 / std::sqrt( 1.0 + norm_2_sq(temp) );
      cs = temp * sn;
   } 
   else 
   {
      Real temp = dy / dx;
      cs = 1.0 / std::sqrt( 1.0 + norm_2_sq(temp) );
      sn = temp * cs;
   }

#if 0
   // debug: This is the matrix that we want to be Unitary
   LinearAlgebra::Matrix<std::complex<double> > M(2,2);
   M(0,0) = cs;
   M(1,1) = cs;
   M(0,1) = conj(sn) * (cs / conj(cs));
   M(1,0) = -sn;

   TRACE(M*herm(M));
#endif
}

template<typename Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, Real cs, Real sn)
{
   Real temp  =  conj(cs) * dx + conj(sn) * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}

template <typename Vector, typename MultiplyFunc, typename PrecFunc>
int 
GmRes(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int& m, int& max_iter, double& tol, PrecFunc Precondition, int Verbose = 0)
{
  //  typedef typename Vector::value_type value_type;
  typedef std::complex<double> value_type;
  typedef LinearAlgebra::Vector<value_type> VecType;
  VecType s(m+1), cs(m+1), sn(m+1);
  LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);

  double normb = norm_frob(Precondition(b));
  Vector r = Precondition(b - MatVecMultiply(x));
  double beta = norm_frob(r);

  // Initial guess for x:
  // let r = MatVecMultiply(x);
  // minimize choose prefactor to minimize norm_frob_sq(b - a*r)
  // This means 
  // norm_frob_sq(b) + |a|^2 norm_frob_sq(r) - <b|ar> - <ar|b>
  // let <b|r> = c
  // then norm_frob_sq(b) + |a|^2 norm_frob_sq(r) - 2*real(ac)
  // So choose a = k*conj(c), for k real
  // Then minimize:
  // f(k) = norm_frob_sq(b) + k^2 |c|^2 norm_frob_sq(r) - 2*k*|c|^2
  // minimize with respect to k:
  // df/dk = 0 =>
  // 2k |c|^2 norm_frob_sq(r) - 2|c|^2 = 0
  // -> k = 1/norm_frob_sq(r)
  // so scale r => r*conj(<b|r>)/norm_frob(r)
  // r -> <r|b>/norm_frob(r)
  // In practice, this seems to have negligble effect versus simply scaling
  // the guess vector by the norm of b, so it is currently disabled.
  // Also I don't understand the complex conjugation here.
  
#if 0
  TRACE(norm_frob(r));
  r = MatVecMultiply(x);
  r = Precondition(b - (inner_prod(b,r)/norm_frob_sq(r))*r);
  beta = norm_frob(r);
  TRACE(norm_frob(r));
#endif

  if (normb == 0.0)
    normb = 1;
  double resid = norm_frob(r) / normb;
  DEBUG_TRACE(norm_frob(b))(norm_frob(MatVecMultiply(x)));
  DEBUG_TRACE(resid)(norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
  if ((resid = norm_frob(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector w;

  Vector* v = new Vector[m+1];

  int j = 1;
  while (j <= max_iter)
  {
     v[0] = (1.0 / beta) * r;
     zero_all(s);
     s[0] = beta;

     int i = 0;
     while (i < m && j <= max_iter && resid >= tol)
     {
	if (Verbose > 2)
	   std::cerr << "GMRES: iteration " << i << std::endl;

	double NormFrobSqH = 0; // for DGKS correction

        w = Precondition(MatVecMultiply(v[i]));
        for (int k = 0; k <= i; k++) 
        {
           H(k, i) = inner_prod(v[k], w);
           w -= H(k, i) * v[k];
	   NormFrobSqH += LinearAlgebra::norm_frob_sq(H(k,i));
        }

	// Apply DGKS correction, if necessary
	double NormFrobSqF = norm_frob_sq(w);
	if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
	{
	   DEBUG_TRACE("DGKS correction in GMRES")(NormFrobSqF / (DGKS_Threshold * DGKS_Threshold * NormFrobSqH));
	   for (int k = 0; k <= i; k++) 
	   {
	      value_type z = inner_prod(v[k], w);
	      H(k, i) += z;
	      w -= z * v[k];
	   }
	}

	// Continue with our normal schedule...
        H(i+1, i) = norm_frob(w);
        v[i+1] = w * (1.0 / H(i+1, i));

        for (int k = 0; k < i; k++)
           ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
        GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
	resid = norm_2(s[i+1]) / normb;

	if (Verbose > 2)
	   std::cerr << "GMRES: resid=" << resid << '\n';

	++i;
	++j;

#if 0
	// debugging check
	{
	   Vector X2 = x;
	   Update(X2, i-1, H, s, v);
	   Vector R = Precondition(b - MatVecMultiply(X2));
	   TRACE(i)(norm_2(s[i]))(norm_frob(R));
	}
#endif

     }
     Update(x, i-1, H, s, v);
     r = Precondition(b - MatVecMultiply(x));
     beta = norm_frob(r);

     // use the old value of resid here, to avoid cases 
     // where the recalculation no longer satisfies the convergence criteria
     if (resid < tol) 
     {
	double UpdatedResid = beta / normb;
        tol = UpdatedResid;
        max_iter = j;
        delete [] v;
	if (Verbose)
	   std::cerr << "GMRES: finished, iter=" << (j-1) << ", approx resid=" << resid 
		     << ", actual resid=" << UpdatedResid << std::endl;
        return 0;
     }
     else
     {
	if (Verbose > 1)
	   std::cerr << "GMRES: restarting, iter=" << (j-1) << ", resid=" << resid << '\n';
     }
     resid = beta / normb;
     DEBUG_TRACE(resid)(norm_frob(Precondition(b - MatVecMultiply(x))) / normb)
        (norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
  }
  
  tol = resid;
  delete [] v;
  return 1;
}

#endif
