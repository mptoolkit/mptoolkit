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

template < typename Matrix, typename Vector1, typename Vector2 >
void 
Update(Vector1 &x, int k, Matrix &h, Vector2 &s, Vector1 v[])
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
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
   //Real temp  =  cs * dx + conj(sn) * dy * (cs / conj(cs));
   //   TRACE(LinearAlgebra::norm_frob_sq(cs/conj(cs)));
   Real temp  =  cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}

template <typename Vector, typename MultiplyFunc, typename PrecFunc>
int 
GmRes(Vector &x, MultiplyFunc MatVecMultiply, Vector const& b,
      int& m, int& max_iter, double& tol, PrecFunc Precondition, int Verbose = 0)
{
   double resid;
  int i, j = 1, k;
  Vector w;

  //  typedef typename Vector::value_type value_type;
  typedef std::complex<double> value_type;
  typedef LinearAlgebra::Vector<value_type> VecType;
  VecType s(m+1), cs(m+1), sn(m+1);
  LinearAlgebra::Matrix<value_type> H(m+1, m+1, 0.0);
  
  double normb = norm_frob(Precondition(b));
  Vector r = Precondition(b - MatVecMultiply(x));
  double beta = norm_frob(r);
  
  if (normb == 0.0)
    normb = 1;
  resid = norm_frob(r) / normb;
  DEBUG_TRACE(resid)(norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
  if ((resid = norm_frob(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector* v = new Vector[m+1];

  while (j <= max_iter) 
  {
     v[0] = (1.0 / beta) * r;
     zero_all(s);
     s[0] = beta;
    
     for (i = 0; i < m && j <= max_iter; i++, j++) 
     {
	if (Verbose > 1)
	   std::cerr << "GMRES iteration " << i << std::endl;

        w = Precondition(MatVecMultiply(v[i]));
        for (k = 0; k <= i; k++) 
        {
           H(k, i) = inner_prod(w, v[k]);
           w -= H(k, i) * v[k];
        }
        H(i+1, i) = norm_frob(w);
        v[i+1] = w * (1.0 / H(i+1, i)); // ??? w / H(i+1, i)

        for (k = 0; k < i; k++)
           ApplyPlaneRotation(H(k,i), H(k+1,i), cs[k], sn[k]);
      
        GeneratePlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        ApplyPlaneRotation(H(i,i), H(i+1,i), cs[i], sn[i]);
        //TRACE(H(i,i))(H(i+1,i));

        ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      
        if ((resid = norm_2(s[i+1]) / normb) < tol) 
        {
           Update(x, i, H, s, v);
           tol = resid;
           max_iter = j;
           delete [] v;
	   if (Verbose)
	      std::cerr << "GMRES finished, iter=" << (j-1) << ", resid=" << resid << std::endl;

	   DEBUG_TRACE("GMRES return")(resid);
           return 0;
        }
        //TRACE(resid);
     }
     //TRACE(H);
     Update(x, i - 1, H, s, v);
     r = Precondition(b - MatVecMultiply(x));
     beta = norm_frob(r);
     if ((resid = beta / normb) < tol) 
     {
        tol = resid;
        max_iter = j;
        delete [] v;
	if (Verbose)
	   std::cerr << "GMRES finished, iter=" << (j-1) << ", resid=" << resid << std::endl;
	DEBUG_TRACE("GMRES return")(resid);
        return 0;
     }
     else
     {
	if (Verbose)
	   std::cerr << "GMRES restarting, iter=" << (j-1) << std::endl;
     }
     DEBUG_TRACE(resid)(norm_frob(Precondition(b - MatVecMultiply(x))) / normb)
        (norm_frob(b - MatVecMultiply(x)) / norm_frob(b));
  }
  
  tol = resid;
  delete [] v;
  return 1;
}

#endif
