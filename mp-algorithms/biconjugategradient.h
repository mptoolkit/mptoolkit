// -*- C++ -*- $Id$
//
//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// CG follows the algorithm described on p. 15 in the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;

template <typename Vector, typename MultiplyFunctor, typename MultiplyHermFunctor,
          typename PreFunctor, typename HermPreFunctor>
void
BiConjugateGradient(Vector &x, MultiplyFunctor MatVecMultiply, 
                    MultiplyHermFunctor HermVecMultiply,
                    Vector const& b,
                    int& max_iter, double& tol,
                    PreFunctor Precondition, 
                    HermPreFunctor HermPrecondition)
{
  double resid;
  Vector p, z, q;
  Vector ptilde, ztilde, qtilde;
  typedef std::complex<double> value_type;
  value_type alpha, beta, rho, rho_1;

  double normb = norm_frob(b);
  Vector rtilde = b - MatVecMultiply(x);
  Vector r = Precondition(conj(rtilde));

  if (normb == 0.0) 
     normb = 1;
  
  if ((resid = norm_frob(rtilde) / normb) <= tol) 
  {
     tol = resid;
     max_iter = 0;
     return;
  }

  TRACE(resid)(norm_frob(b - MatVecMultiply(x)) / normb);

  for (int i = 1; i <= max_iter; i++) 
  {
     z = Precondition(rtilde);
     ztilde = HermPrecondition(r);

     rho = inner_prod(conj(z), r);
     rho = conj(rho);
    
     if (i == 1)
     {
        p = z;
        ptilde = ztilde;
     }
     else 
     {
        beta = rho / rho_1;

        beta = conj(beta);

        p = z + beta * p;
        ptilde = ztilde + beta * ptilde;
     }
    
     qtilde = MatVecMultiply(p);
     q = HermVecMultiply(ptilde);
     alpha = rho / inner_prod(conj(ptilde), qtilde);
    
     alpha = conj(alpha);

     x += alpha * p;
     r -= alpha * q;
     rtilde -= alpha * qtilde;
    
     if ((resid = norm_frob(rtilde) / normb) <= tol) 
     {
        tol = resid;
        max_iter = i;
      return;     
    }

     TRACE(resid)(norm_frob(b - MatVecMultiply(x)) / normb);

    rho_1 = rho;
  }
  
  tol = resid;
  // at this point, we did not converge
  return;
}
